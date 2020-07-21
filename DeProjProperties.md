# DeProj output and properties.

We give here the list of properties stored in the DeProj classes and their definition.

[TOC]

## The core classes: `deproj` and `epicell`.

The DeProj toolbox revolves around two MATLAB classes:

- The `epicell`class is used to store the data for one cell in a tissue. It is made of several fields we describe below:

```matlab
o = 

  epicell with properties:

                 boundary: [26×3 single]
                   center: [2.4705 11.1826 3.1008]
             junction_ids: [5×1 double]
                     area: 8.0176
                perimeter: 12.5227
             euler_angles: [-2.0734 0.4195 -0.2500]
               curvatures: [0.0110 -4.8103e-05 0.0240 -0.0020]
              ellipse_fit: [2.2284 11.1114 3.1008 2.3848 1.1178 0.4528]
             eccentricity: 0.8834
           proj_direction: 1.2539
         uncorrected_area: 7.3173
    uncorrected_perimeter: 12.1116
                       id: 4
```

- The `deproj` class manages a collection of `epicell`s, and represent for instance the results of the analysis of a whole image.

```matlab
dpr = 

  deproj with properties:

          epicells: [426×1 epicell]
    junction_graph: [1×1 graph]
             units: 'µm'
```

## The `epicell` class properties.

### `boundary`

The `boundary` property stores a 3D polygon that delimits the cell. It is represented as a `N x 3` matrix, with points sorted along the polygon. The coordinates are in physical units. The Z-values are the position of the vertices on the tissue surface, given the by the smoothed height-map. 

There are fewer than on vertex per pixel on the original segmentation image: we prune redundant points to lower the memory footprint of segmentation results.

```matlab
o = dpr.epicells(1);
p = o.boundary;
plot3( p(:,1), p(:,2), p(:,3), 'k-o' )
xlabel('x (µm)'), ylabel('y (µm)'), zlabel('z (µm)')
axis equal
```

<img src="static/BoundaryPlot.png" alt="BoundaryPlot"  width="400" />

### `center`

The cell center in physical coordinates. `1 x 3` array.

### `junction_ids`

The list of indices of junctions that this cell touches.

You can notice that a `deproj` instance has a `junction_graph` property. This is the graph of the junctions in the tissues. A junction is a location on the tissue surface where at least 3 cells connect. We give more details in the next section on `deproj` properties.

```matlab
dpr = 

  deproj with properties:

          epicells: [426×1 epicell]
    junction_graph: [1×1 graph]
             units: 'µm'

>> ids = o.junction_ids

ids =

     1
     5
    10
    16

>> j = dpr.junction_graph.Nodes( ids, : ) 

j =

  4×2 table

             Centroid              ID
    ___________________________    __

    0.366     23.424    0.36189     1
    0.976     24.034     0.3217     5
    1.647     21.228    0.90516    10
    2.257     21.838    0.92981    16

>> pj = j.Centroid;
>> plot3( p(:,1), p(:,2), p(:,3), 'k-o' )
>> axis equal, hold on
>> plot3( pj(:,1), pj(:,2), pj(:,3), 'ro', 'MarkerFaceColor', 'r' )
>> xlabel('x (µm)'), ylabel('y (µm)'), zlabel('z (µm)')
```

<img src="static/BoundaryJunctionPlug.png" alt="BoundaryJunctionPlug" width="400" />

### `area`

The area of the apical surface of the cell, delimeted by its boundary. This value reports the 3D area, of the de-projected cell on the tissue surface. Units are physical units; if you used a pixel size in µm, the area will be in µm².

### `perimeter`

The perimeter of the cell boundary, in physical coordinates. Again, using the real 3D coordinates and physical units of length.

### `euler_angles`

The orientation of the cell apical surface plane. This property is not entirely trivial.

The 3D boundary of the cell define an oblique plane to which they are the closest. This local plane has a certain orientation, that we give as the Euler angles, in radians, as a `1 x 3` array. The Euler angles are reported using the [ZX'Z'' convention](https://en.wikipedia.org/wiki/Euler_angles#Chained_rotations_equivalence). Repeating what we said in the example:

- The first one, `alpha` is the orientation of the cell plane. As an analogy, imaging you are facing a hill, the slope going up. The direction (south, west…) in which the slope is the largest is given by the angle `alpha`. Values within `[-π ; π[`. 
- The second one, `beta` measures the slope of this plane with XY plane (middle panel). A value of 0º indicates that the cell plane is parallel to XY. Values within `[0 ; π[`. 
- The third one , `gamma` measures the cell main orientation in the cell plane (bottom panel). Because the cell plane was rotated a first time by `alpha`, this angle does not give a result immediately usable. 

### `curvatures`

Local curvature metrics at the cell center. 

The height-maps gives the shape of the tissue surface, so we can derive the local curvature from it. The `curvatures` property is a `1 x 4` array that gives respectively:

- The [mean curvature](https://en.wikipedia.org/wiki/Mean_curvature).
- The [Gaussian curvature](https://en.wikipedia.org/wiki/Gaussian_curvature).
- The first [principle curvature](https://en.wikipedia.org/wiki/Principal_curvature).
- The second principle curvature. 

They are in physical units. If you used a pixel size in µm, the mean curvature, the first and second principle cuvature will be in 1/µm, and the Gaussian curvature will be in 1/µm².

### `ellipse_fit`

Results of the fit of a 2D ellipse on the 3D boundary of the cell.

The cell apical plane we discribe above, and which orientation is given by the `euler_angles` property, is the plane to which the 3D points of the boundary are the closest. If we project the 3D points on this plane, we can fit an ellipse to the projected 2D plane and yield a description of the cell extend and orientation. Caution: the ellipse fit is not made on the XY plane, but on an oblique plane locally tangent to the tissue. Check the Figure 3 on the example, and rotate it, to see that the ellipses have a 3D orientation.

The `ellipse_fit` contains the result of the fit, in this tangent plane, as a `1 x 6` array, containing respectively:

- The x coordinate of the ellipse center, in physical units.
- The y coordinate of the ellipse center, in physical units.
- The z coordinate of the ellipse center, in physical units. These 3 coordinates should have values very close that of the `center` property.
- The [semi-major axis of the ellipse](https://en.wikipedia.org/wiki/Ellipse#Semi-major_and_semi-minor_axes) `a`, in physical units.
- The semi-minor axis of the ellipse `b`, in physical units. We always have `a > b` and both values are positive.
- The angle of the semi-major axis with the X'' axis in the rotated tangent plane, in radians. Because this plane was rotated twice with respect to the tissue (euler angles alpha and gamma, this value is not directly usable. See the `proj_direction` property.

### `eccentricity`

The ellipse [eccentricity](https://en.wikipedia.org/wiki/Ellipse#Eccentricity), derived from the `ellipse_fit` field. 

It measures how elongated it the cell and varies from 0 to 1. Cells with an eccentricity of 0 resemble a circle. Cells with an eccentricity close to 1 are very elongated.

### `proj_direction`

The angle of the long axis of the cell ellipse, with respect to the X axis of the tissue. In radians. See the Figure 3 of the main example.

### `uncorrected_area`

The area of the cell, *if it was projected and measured on the XY plane*.This property is included only to assess the impact of the projection distorsion artifact. In physical units.

### `uncorrected_perimeter`

The same, for the cell perimeter.

### `id`

The cell unique ID within a `deproj` instance. Strictly positive integers.

## The `deproj` class properties.

A `deproj` instance has only 3 properties.

### `epicells`

A `N x 1` array of `epicell` instances. Each `epicell` represents a cell in the tissue. See above for its properties. The `id` property of an `epicell` corresponds to its index in this array.

### `junction_graph`

A undirected graph of junction connections.

 A junction is a location on the tissue surface where at least 3 cells connect.

 <img src="static/JunctionExample.png" alt="JunctionExample" width="300" />

The graph stored in the `junction_graph` property is made of the junctions as nodes, and of their connection of edges. If two junctions are linked by an edge in this graph, it means that there is exactly one ridge between two cells that connects them.

The nodes store their ID and the 3D position of their centroid, in physical units. Because a junction can be made of several pixels, the `Centroid` matrix gives the X, Y, Z position of the pixels that compose a junction.

```matlab
>> g = dpr.junction_graph

g = 

  graph with properties:

    Edges: [1329×1 table]
    Nodes: [920×2 table]

>> head(g.Nodes)

ans =

  8×2 table

             Centroid              ID
    ___________________________    __

    0.366     23.424    0.36189    1 
    0.732     19.581    0.98203    2 
    0.732     24.705    0.19083    3 
    1.098      4.575      3.898    4 
    0.976     24.034     0.3217    5 
    1.281     10.614     2.8064    6 
    1.464      3.294     4.3027    7 
    1.586     14.335     2.0956    8 
```

The edges simply store theirs source and target node ids:

```matlab
>> head(g.Edges)

ans =

  8×1 table

    EndNodes
    ________

    1     5 
    1    10 
    3     5 
    3    11 
    4     7 
    4    14 
    5    16 
    6    12 
```

Edges are undirected and there are no duplicates.

We can use this junction graph to plot a topological representation of the segmentation results:

```matlab
% Path to the images.
>> root_folder = 'samples';

% Load the segmentation image.
>> mask_filename       = 'Segmentation-2.tif';
>> I = imread( fullfile( root_folder, mask_filename ) );

% Because we want to display the image with the same coordinates
% that of the junctions, we need to specify the pixel size and use
% the 'XData' and 'YData' arguments of imshow.
>> pixel_size = 0.183; % µm

% Display the segmentation.
>> imshow( ~I , [ 0 1 ], ...
	'Border', 'tight',  ...
	'Xdata', pixel_size * [ 1 size(I,2) ], ...
	'YData', pixel_size * [ 1 size(I,1) ] )
>> hold on

% Get the junction graph (run the RunExample.m file first)
g = dpr.junction_graph;

% Display the junction graph. We use the graph builtin plot method.
>> plot( g, ...
    'XData', g.Nodes.Centroid(:,1), ...
    'YData', g.Nodes.Centroid(:,2), ...
    'LineWidth', 2, ...
    'EdgeColor', 'b', ...
    'EdgeAlpha', 1, ...
    'Marker', 'o', ...
    'MarkerSize', 4, ...
    'NodeColor', 'r' )
```

![PlotJunctionGraph_1](static/PlotJunctionGraph_1.png)

![PlotJunctionGraph_2](static/PlotJunctionGraph_2.png)

### `units`

This property simply stores the name of the physical unit of length that was specified at creation. It is useful to generate properly annotated figures.

