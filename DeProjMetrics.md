# DeProj metrics.

We give here the list of metrics returned by DeProj and their definition.

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

The `boundary` property stores a 3D polygon that delimits the cell. It is represented as a `N x 3` matrix, with points sorted 

## The `deproj` class properties.