[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build status](https://github.com/PetrKryslUCSD/ShellStructureTopo.jl/workflows/CI/badge.svg)](https://github.com/PetrKryslUCSD/ShellStructureTopo.jl/actions)

# ShellStructureTopo.jl

The purpose of this package is to detect topological classification of a
triangular surface mesh. For instance, the orientation of the triangles can be
unified across the surface mesh. The topological faces can be identified.

![Sample of unified normals](sample_normals50.png)
![Sample of detected topological faces](sample_topology50.png)

## Usage

### When used with MeshCore

Given an incidence relation, `t2v` (triangles to vertices),
we can produce an oriented surface: 
```
orientedt2v, orientable = orient_surface_mesh(t2v)
```
And, when the triangles tile a single surface, we can check that all the triangles have been classified on surface 1:
```
@test length(unique(attribute(orientedt2v.left, "surfid"))) == 1
```


Given an incidence relation, `t2v` (triangles to vertices),
we can classify the triangles on oriented surfaces: 
```
t2v = make_topo_faces(t2v)
```
In this instance, we would expect the triangles to represent four distinct topological surfaces. We can check the `"tfid"` attribute of the incidence relation `t2v`:
```
@test length(unique(attribute(t2v.left, "tfid"))) == 4
```

The topological classification may be visualized with
```
t2v = make_topo_faces(t2v)
vtkwrite("mt008-classified", t2v, [(name = "tfid", )])
```
The topological surfaces will be labeled with the attribute `"tfid"`.

### When used with FinEtools

The boundary triangles of the triangulated surface of a solid part may be classified into topological faces with:
```
fens, bfes = make_topo_faces(fens, bfes)
```
Here `bfes.label` records the numbers of the topological surfaces.

The classification may be visualized with
```
VTK.vtkexportmesh("mt013.vtk", connasarray(bfes), fens.xyz, VTK.T3; scalars=[("topological_face", bfes.label)]);
```

The library may be also used to create partition of the topological faces into individual patches.  
```
surfids, partitionids = create_partitions(fens, fes, 50)
```

The classification into topological faces and the partitioning may then be viewed with
```
VTK.vtkexportmesh("mt016_part.vtk", connasarray(fes), fens.xyz, VTK.T3; scalars=[("topological_face", surfids), ("partitioning", partitionids)]);
```

## News

- 07/06/2024: Improve estimation of number of partitions.
- 07/06/2024: Update for FinEtools 8.
