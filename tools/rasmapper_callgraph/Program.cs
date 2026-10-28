// Every method reachable from a RASMapper method, read from the IL of HEC's
// DLLs (calls, constructors, delegates; virtual calls resolved to types that
// are created along the way or passed as already existing).
//
//   CALLGRAPH_STOP=stop.tsv dotnet run -c Release -- <HEC-RAS dir> \
//       RasMapperLib.RASD2FlowArea CreatePropertyTables reach.tsv 4 \
//       RasMapperLib.TerrainLayer RasMapperLib.LandCoverLayer ...
//
// Writes reach.tsv (methods followed) and reach.tsv.stopped (methods the
// rules in stop.tsv record but don't follow: UI, progress, I/O, project model,
// other features).
using Mono.Cecil;
using Mono.Cecil.Cil;

var dir = args[0];
var rootType = args[1];
var rootMethod = args[2];
var outPath = args[3];

string[] hec = { "RasMapperLib", "Geospatial.Core", "Geospatial.GDALAssist", "Geospatial.IO",
                 "Utility.Core", "H5Assist", "TiffAssist", "TiffBinaryReader", "clipper_library" };
var resolver = new DefaultAssemblyResolver();
resolver.AddSearchDirectory(dir);
var rp = new ReaderParameters { AssemblyResolver = resolver };
var asms = hec.Select(n => Path.Combine(dir, n + ".dll")).Where(File.Exists)
              .Select(p => AssemblyDefinition.ReadAssembly(p, rp)).ToList();
var hecNames = asms.Select(a => a.Name.Name).ToHashSet();

IEnumerable<TypeDefinition> AllTypes(TypeDefinition t) =>
    new[] { t }.Concat(t.NestedTypes.SelectMany(AllTypes));
var types = asms.SelectMany(a => a.MainModule.Types).SelectMany(AllTypes).ToList();
var methods = types.SelectMany(t => t.Methods).ToList();
var byName = methods.Where(x => x.HasBody)
    .GroupBy(x => x.Name.Split('.').Last() + "/" + x.Parameters.Count)
    .ToDictionary(g => g.Key, g => g.ToList());
var typeByName = new Dictionary<string, TypeDefinition>();
foreach (var t in types) typeByName.TryAdd(t.FullName, t);
var derivesCache = new Dictionary<(string, string), bool>();

// implementations of a virtual / interface method among the loaded types
var implCache = new Dictionary<string, List<MethodDefinition>>();
List<MethodDefinition> Implementations(MethodDefinition m)
{
    var key = m.FullName;
    if (implCache.TryGetValue(key, out var hit)) return hit;
    var res = new List<MethodDefinition>();
    if (byName.TryGetValue(m.Name.Split('.').Last() + "/" + m.Parameters.Count, out var cands))
        foreach (var cand in cands)
            if (IsOrDerives(cand.DeclaringType, m.DeclaringType.FullName) || cand.Overrides.Any(o => o.FullName == m.FullName))
                res.Add(cand);
    implCache[key] = res;
    return res;
}

var nParams = args.Length > 4 ? int.Parse(args[4]) : -1;
var root = types.First(t => t.FullName == rootType).Methods
                .Where(m => m.Name == rootMethod && m.HasBody && (nParams < 0 || m.Parameters.Count == nParams)).ToList();
// Types that exist before the table step runs (the project's layers); their
// virtual methods can be the target of calls. Extra seeds: args[5..]
var instantiated = new HashSet<string>(args.Skip(5));
var virtualSites = new List<(MethodDefinition md, string from, int depth)>();
bool IsOrDerives(TypeDefinition t, string baseName)
{
    var key = (t.FullName, baseName);
    if (derivesCache.TryGetValue(key, out var hit)) return hit;
    bool found = false;
    for (var b = t; b != null && !found; )
    {
        if (b.FullName == baseName) found = true;
        else if (b.Interfaces.Any(i => { try { return i.InterfaceType.Resolve()?.FullName == baseName; } catch { return false; } })) found = true;
        else { try { b = b.BaseType?.Resolve(); } catch { b = null; } }
    }
    derivesCache[key] = found;
    return found;
}
bool InstantiatedSubtypeHas(MethodDefinition impl)
{
    foreach (var name in instantiated)
        if (typeByName.TryGetValue(name, out var t) && IsOrDerives(t, impl.DeclaringType.FullName)) return true;
    return false;
}
var seen = new Dictionary<string, (MethodDefinition m, int depth, string from)>();
var queue = new Queue<(MethodDefinition m, int depth, string from)>();
foreach (var r in root) { seen[r.FullName] = (r, 0, "(root)"); queue.Enqueue((r, 0, "(root)")); }

// Subsystems that can't change table values: reached, recorded, not followed.
var stopFile = Environment.GetEnvironmentVariable("CALLGRAPH_STOP");
var stops = stopFile != null ? File.ReadAllLines(stopFile).Where(l => l.Trim().Length > 0 && !l.StartsWith("#"))
                                  .Select(l => l.Split('\t')).Select(a => (cat: a[0], rx: new System.Text.RegularExpressions.Regex(a[1]))).ToList()
                             : new List<(string cat, System.Text.RegularExpressions.Regex rx)>();
var stopped = new Dictionary<string, (MethodDefinition m, string cat, string from)>();
string? StopCategory(MethodDefinition m)
{
    var name = m.DeclaringType.FullName + "::" + m.Name;
    foreach (var (cat, rx) in stops) if (rx.IsMatch(name)) return cat;
    return null;
}

void Visit(MethodDefinition target, int depth, string from)
{
    if (target == null || !hecNames.Contains(target.Module.Assembly.Name.Name)) return;
    if (seen.ContainsKey(target.FullName) || stopped.ContainsKey(target.FullName)) return;
    var cat = StopCategory(target);
    if (cat != null) { stopped[target.FullName] = (target, cat, from); return; }
    seen[target.FullName] = (target, depth, from);
    queue.Enqueue((target, depth, from));
    var cctor = target.DeclaringType.Methods.FirstOrDefault(x => x.IsConstructor && x.IsStatic);
    if (cctor != null) Visit(cctor, depth, target.FullName);
}

while (queue.Count > 0)
{
    var (m, depth, _) = queue.Dequeue();
    if (!m.HasBody) continue;
    foreach (var ins in m.Body.Instructions)
    {
        if (ins.Operand is not MethodReference mr) continue;
        MethodDefinition? md = null;
        try { md = mr.Resolve(); } catch { }
        if (md == null) continue;
        if (ins.OpCode == OpCodes.Newobj && hecNames.Contains(md.Module.Assembly.Name.Name)
            && instantiated.Add(md.DeclaringType.FullName))
        {
            // a newly created type may be the target of virtual calls already seen
            foreach (var (vmd, vfrom, vdepth) in virtualSites.ToList())
                foreach (var impl in Implementations(vmd).Where(i => IsOrDerives(md.DeclaringType, i.DeclaringType.FullName)))
                    Visit(impl, vdepth, vfrom + " [virtual]");
        }
        if ((ins.OpCode == OpCodes.Callvirt || ins.OpCode == OpCodes.Ldvirtftn) && (md.IsVirtual || md.IsAbstract))
        {
            virtualSites.Add((md, m.FullName, depth + 1));
            if (!md.IsAbstract) Visit(md, depth + 1, m.FullName);
            foreach (var impl in Implementations(md).Where(InstantiatedSubtypeHas))
                Visit(impl, depth + 1, m.FullName + " [virtual]");
        }
        else Visit(md, depth + 1, m.FullName);
    }
}

using var w = new StreamWriter(outPath);
w.WriteLine("depth\tassembly\ttype\tmethod\til_bytes\treached_from");
foreach (var (m, depth, from) in seen.Values.OrderBy(v => v.depth).ThenBy(v => v.m.DeclaringType.FullName).ThenBy(v => v.m.Name))
    w.WriteLine($"{depth}\t{m.Module.Assembly.Name.Name}\t{m.DeclaringType.FullName}\t{m.Name}({string.Join(",", m.Parameters.Select(p => p.ParameterType.Name))})\t{(m.HasBody ? m.Body.CodeSize : 0)}\t{from}");
using (var sw = new StreamWriter(outPath + ".stopped"))
{
    sw.WriteLine("category\tassembly\ttype\tmethod\til_bytes\treached_from");
    foreach (var (m, cat, from) in stopped.Values.OrderBy(v => v.cat).ThenBy(v => v.m.DeclaringType.FullName))
        sw.WriteLine($"{cat}\t{m.Module.Assembly.Name.Name}\t{m.DeclaringType.FullName}\t{m.Name}\t{(m.HasBody ? m.Body.CodeSize : 0)}\t{from}");
}
foreach (var g in stopped.Values.GroupBy(v => v.cat).OrderBy(g => g.Key))
    Console.WriteLine($"  stopped at {g.Key}: {g.Count()} entry methods");
Console.WriteLine($"{seen.Count} methods reachable from {rootType}.{rootMethod}");
foreach (var g in seen.Values.GroupBy(v => v.m.Module.Assembly.Name.Name).OrderByDescending(g => g.Count()))
    Console.WriteLine($"  {g.Key}: {g.Count()} methods, {g.Sum(v => v.m.HasBody ? v.m.Body.CodeSize : 0)} IL bytes");
