## Synopsis

This program aims to provide an intuitive visual explanation of **Principal Component Analysis** or **PCA**.

Principal Component Analysis is generally used to **reduce the dimension** of a dataset. In a sense, it is a 'compression' of the data. 

If our data is 3D data- for example spatial data -then we can explore a PCA visually. This is because reducing the dimension of 3D data is equivalent to 'flattening it'.

One way of 'flattening' data is simply to ignore or discard one dimension of the data. PCA, however, tells us a way to orientate our data such that we can discard a dimension of the data without losing too much information. The best way to understand this idea is to run this program yourself.

### Screen Recording

https://vimeo.com/156145588

### Bounding Ellipsoid Construction

<a href="https://vimeo.com/156145588">
    <img src="./ellipsoid.png" alt="ellipsoid" width="500px"/>
</a>

### Dimensionality Reduction

<a href="https://vimeo.com/156145588">
    <img src="./dim_reduced.png" alt="dimensionality_reduction" width="500px"/>
</a>

<br><br>

## Motivation

* To explore visual communication of technical ideas in maths and data science
* To create a solid Matrix.class for future computer graphics and data analysis projects
* To learn and explore Linear Algebra 

## Installation

* Download and install Processing from http://www.processing.org.
* Ensure all files are contained within a directory called **Principal_Component_Analysis**

## Related Work

* Visualisation of the underlying algorithm for PCA:

https://github.com/cbrookhouse8/Jacobi_Rotations_Visualised

## Futher development

* Implement vector2, vector3, vector4 classes rather that always using the matrix class.
* Explore recursive matrix multiplication
* Refactor Quaternion and Matrix classes as classes within a Java library that is imported into Processing.
* Rename 'Project' button. Or at least specify what kind of 'projection' this is
* Improve camera mouse dragging
* Create a dedicated class for the Camera, with associated orientation properties etc
* Create a small set of parsable text instructions that can coordinate the animation and transitions
* Add text to indicate exactly what is happening in each stage of PCA

## References

Eric Lengyel, *Mathematics for 3D Game Programming and Computer Graphics*

Mark Richardson, *Principal Component Analysis*

