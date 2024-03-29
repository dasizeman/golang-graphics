package graphics

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"log"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"

	"github.com/gonum/matrix/mat64"
)

/* ----Structures---- */

/* RGBColor */

// RGBColor is an RGB color with 8bit channels
type RGBColor struct {
	R, G, B byte
}

/* SoftFrameBuffer */

// SoftFrameBuffer is a simulated frame buffer that may be written to a file
type SoftFrameBuffer struct {
	Width, Height int
	buffer        [][]RGBColor
	zBuffer       [][]float64
	colors        map[RGBColor]string
}

// NewSoftFrameBuffer allocates a new SoftFrameBuffer with specified width and
// height, with 8bit color channels
func NewSoftFrameBuffer(width, height int) *SoftFrameBuffer {
	newBuffer := &SoftFrameBuffer{
		Width:  width,
		Height: height}

	black := RGBColor{0, 0, 0}

	newBuffer.colors = make(map[RGBColor]string)

	// Allocate the pixel buffer
	newBuffer.buffer = make([][]RGBColor, newBuffer.Height)
	newBuffer.zBuffer = make([][]float64, newBuffer.Height)
	for y := range newBuffer.buffer {
		newBuffer.buffer[y] = make([]RGBColor, newBuffer.Width)
		newBuffer.zBuffer[y] = make([]float64, newBuffer.Width)
	}

	// Set the buffer to black and initialize z buffer
	for y := 0; y < newBuffer.Height; y++ {
		for x := 0; x < newBuffer.Width; x++ {
			newBuffer.WritePixel(x, y, black)
			newBuffer.zBuffer[y][x] = -1
		}
	}

	return newBuffer
}

func (frameBuffer *SoftFrameBuffer) writeDepthCueColor(z, front, back float64, x, y int, baseColor RGBColor) {

	//fmt.Fprintf(os.Stderr, "%f\n", z)
	// Calculate the scale factor based on the distance from the back
	scale := (z - back) / (front - back)

	// Scale the base color
	if x == 400 {
		//fmt.Println()
	}
	scaledColor := scaleColor(baseColor, scale)
	if y == 233 && scaledColor.R == 7 {
		//fmt.Printf("%d,%d,%d\n", scaledColor.R, scaledColor.G, scaledColor.B)
		//fmt.Printf("%f\n", z)
	}

	// Write the pixel
	frameBuffer.WritePixel(x, y, scaledColor)

}

func get20ColorXPMCode(color RGBColor) string {
	code := ""
	currentColor := &color.R
	for i := 0; i < 3; i++ {
		// Map the channel value to an ascii value between 48 and 68
		ascii := 48 + (float64(*currentColor)/255)*20
		code += string(int(ascii))

		if currentColor == &color.R {
			currentColor = &color.G
		} else {
			currentColor = &color.B
		}
	}

	return code
}

func scaleColor(color RGBColor, scale float64) RGBColor {
	return RGBColor{
		byte(float64(color.R) * scale),
		byte(float64(color.G) * scale),
		byte(float64(color.B) * scale)}
}

func (frameBuffer *SoftFrameBuffer) updateZBuffer(x, y int, z float64) bool {
	if z <= frameBuffer.zBuffer[y][x] {
		return false
	}
	frameBuffer.zBuffer[y][x] = z
	return true
}

// WritePixel writes a pixel color to the buffer, assuming (0,0) is the
//bottom-left corner
func (frameBuffer *SoftFrameBuffer) WritePixel(x, y int, color RGBColor) {

	// We will reverse the y origin for compatibility with XPM output
	frameBuffer.buffer[frameBuffer.Height-y-1][x] = color

	// Add the color to the buffer's color map
	frameBuffer.colors[color] = get20ColorXPMCode(color)

}

/* Port */

// Port2D is a rectangular geometry.  This can be used as, for example, a world
// window or a viewport
type Port2D struct {
	XMin, YMin, XMax, YMax float64
}

/* Scene */

// Scene is a virtual space that can be rendered to a frame buffer
type Scene struct {
	// Objects in the scene
	Objects []Geometry

	// All vertices
	Vertices []*Vertex

	/*
		// Global scene scale about world origin
		Scale float64

		// Global scene rotation, in degrees, counter-clockwise about origin
		Rotation int

		// Global translation in the X direction
		XTranslation float64

		// Global translation in the Y direction
		YTranslation float64
	*/

	// Projection Reference Point
	PRP *mat64.Vector

	// View Reference Point
	VRP *mat64.Vector

	// View Plane Normal
	VPN *mat64.Vector

	// View Up Vector
	VUP *mat64.Vector

	// View Reference Coordinate window
	ViewPlane Port2D

	// Clipping Planes
	FrontClippingPlane, BackClippingPlane float64

	// Viewport
	Viewport Port2D

	// Projection method
	UseOrthographicProjection bool
}

// Render renders the scene to the given buffer
func (scene *Scene) Render(buffer *SoftFrameBuffer) {
	// Get the view matrix
	viewMatrix := GetViewMatrix(scene.VPN, scene.VUP, scene.VRP)

	var projectionMatrix mat64.Matrix
	if scene.UseOrthographicProjection {
		projectionMatrix = GetOrthographicProjectionMatrix(scene.PRP, scene.ViewPlane,
			scene.FrontClippingPlane, scene.BackClippingPlane)
	} else {
		projectionMatrix = GetPerspectiveProjectionMatrix(scene.PRP, scene.ViewPlane,
			scene.BackClippingPlane)
	}
	//fmt.Fprintf(os.Stderr, "%v\n", mat64.Formatted(viewMatrix))
	//fmt.Fprintf(os.Stderr, "%v\n", mat64.Formatted(projectionMatrix))

	finalMatrix := mat64.NewDense(4, 4, nil)
	finalMatrix.Mul(projectionMatrix, viewMatrix)

	// Transform the viewport for all vertices
	for _, vertex := range scene.Vertices {
		// Create a column vector from the vertex
		vector := mat64.NewVector(len(vertex.attributes), vertex.attributes)

		// Transform the view volume
		vector.MulVec(finalMatrix, vector)

		vertex.attributes = stripVector(vector)
	}

	var filteredPolygons []Geometry
	for _, polygon := range scene.Objects {
		if isPolygonValid(polygon, scene) {
			filteredPolygons = append(filteredPolygons, polygon)
		}
	}

	// Project all the vertices
	for _, vertex := range scene.Vertices {
		// Create a column vector from the vertex
		vector := mat64.NewVector(len(vertex.attributes), vertex.attributes)

		// Do projection if needed
		if !scene.UseOrthographicProjection {
			// For perspective projection, divide x and y by -z
			for dim := 0; dim < 2; dim++ {
				vector.SetVec(dim, vector.At(dim, 0)/-vector.At(2, 0))
			}
		}

		// Put the new vertex data back
		vertex.attributes = stripVector(vector)

		// "Convert" the vertex back to 2D
		//vertex.attributes = vertex.attributes[:3]
		//vertex.attributes[2] = 1

		// Apply world-to-viewport transformation
		canonicalViewport := Port2D{
			XMin: -1,
			YMin: -1,
			XMax: 1,
			YMax: 1}
		scene.transformToViewport(vertex, canonicalViewport)
	}

	// Draw the projected vertices
	for i := range filteredPolygons {

		filteredPolygons[i].Discretize()
		filteredPolygons[i].scanFill(buffer, scene)
		//filteredPolygons[i].Draw(buffer)

	}
}

func isPolygonValid(geo Geometry, scene *Scene) bool {
	var checker func(*Vertex, *Scene) int
	if scene.UseOrthographicProjection {
		checker = getOrthogonalBitCode
	} else {
		checker = getPerspectiveBitCode
	}
	valid := true
	for _, vertex := range geo.vertices {
		if checker(vertex, scene) != 0 {
			valid = false
		}
	}
	return valid
}

func (scene *Scene) transformToViewport(vertex *Vertex, from Port2D) {

	// Get the matrix to translate the world window to the origin
	originTranslationMatrix :=
		GetTranslationTransformMatrix3D(-float64(from.XMin),
			-float64(from.YMin), 0)
	//fmt.Printf("%v\n", mat64.Formatted(originTranslationMatrix))

	// Get the scale matrix
	xScale := float64(scene.Viewport.XMax-scene.Viewport.XMin) /
		float64(from.XMax-from.XMin)

	yScale := float64(scene.Viewport.YMax-scene.Viewport.YMin) /
		float64(from.YMax-from.YMin)

	viewportScaleMatrix :=
		GetScaleTransformMatrix3D(xScale, yScale, 1)
	//fmt.Printf("%v\n", mat64.Formatted(viewportScaleMatrix))

	// Get the matrix to translate the viewport to the desired location
	viewportTranslationMatrix :=
		GetTranslationTransformMatrix3D(float64(scene.Viewport.XMin),
			float64(scene.Viewport.YMin), 0)
	//fmt.Printf("%v\n", mat64.Formatted(viewportTranslationMatrix))

	// Compose the final transformation matrix
	transformationMatrix := mat64.NewDense(4, 4, nil)
	transformationMatrix.Mul(viewportScaleMatrix, originTranslationMatrix)
	transformationMatrix.Mul(viewportTranslationMatrix, transformationMatrix)

	// Apply the transformation
	vertex.Transform(transformationMatrix)

}

/* Vertex */

// Vertex is a vertex with attributes
type Vertex struct {
	attributes  []float64
	Transformed bool
}

// Create2DVertex creates a vertex with the specified x and y coordinates
func Create2DVertex(a, b float64) *Vertex {
	// TODO use this where we create 2D vertices
	polygons := &Vertex{}
	polygons.AddAttribute(a)
	polygons.AddAttribute(b)

	// Homogenous
	polygons.AddAttribute(1)
	return polygons
}

// Create2DVertexInt is a wrapper for Create2DVertex
// that accepts integer arguments
func Create2DVertexInt(a, b int) *Vertex {
	return Create2DVertex(float64(a), float64(b))
}

// Create3DVertex creates a vertex with the specified x,y, and z coordinates
func Create3DVertex(a, b, c float64) *Vertex {
	polygons := &Vertex{}
	polygons.AddAttribute(a)
	polygons.AddAttribute(b)
	polygons.AddAttribute(c)

	// Homogenous
	polygons.AddAttribute(1)

	return polygons
}

// AddAttribute adds an attribute to vertex
func (vertex *Vertex) AddAttribute(attr float64) {
	vertex.attributes = append(vertex.attributes, attr)
}

// AddIntAttribute adds an integer attribute to the vertex
func (vertex *Vertex) AddIntAttribute(attr int) {
	vertex.AddAttribute(float64(attr))
}

func (vertex *Vertex) String() string {
	return fmt.Sprintf("%v", vertex.attributes)
}

// Equal checks if this Vertex is equal to another Vertex
func (vertex *Vertex) Equal(other *Vertex) bool {
	if other == nil {
		return false
	}
	if len(vertex.attributes) != len(other.attributes) {
		return false
	}

	for i := 0; i < len(vertex.attributes); i++ {
		if vertex.attributes[i] != other.attributes[i] {
			return false
		}
	}

	return true
}

// All of these vector operations expect two Vertices of
// the same dimension

// Magnitude returns the vector magnitude of this vertex
func (vertex *Vertex) Magnitude() float64 {
	squareSum := 0.0
	for _, attribute := range vertex.attributes {
		squareSum += math.Pow(attribute, 2)
	}

	return math.Sqrt(squareSum)
}

// Normalized returns the normalized vector of this vertex
func (vertex *Vertex) Normalized() *Vertex {
	res := &Vertex{}
	magnitude := vertex.Magnitude()
	for _, attribute := range vertex.attributes {
		res.AddAttribute(attribute / magnitude)
	}

	return res
}

// TwoDNormal returns the normalized vector normal to the 2D vertex
func (vertex *Vertex) TwoDNormal() *Vertex {
	normalized := vertex.Normalized()
	res := &Vertex{}

	// x = y
	// y = -x
	res.AddAttribute(normalized.attributes[1])
	res.AddAttribute(-1 * normalized.attributes[0])

	return res
}

// VAdd vector adds two Vertices
func VAdd(a, b *Vertex) *Vertex {
	res := &Vertex{}
	for i := 0; i < len(a.attributes); i++ {
		res.AddAttribute(a.attributes[i] + b.attributes[i])
	}

	return res
}

// VScale Vector scales vertex by factor
func VScale(vertex *Vertex, factor float64) *Vertex {
	res := &Vertex{}
	for _, attribute := range vertex.attributes {
		res.AddAttribute(attribute * factor)
	}

	return res
}

// VDot returns the vector dot product of the vertices
func VDot(a, b *Vertex) float64 {
	sum := 0.0
	for i := 0; i < len(a.attributes); i++ {
		sum += a.attributes[i] * b.attributes[i]
	}

	return sum
}

/* Tranformations */

// GetScaleTransformMatrix3D returns a transformation matrix for a scale
// transform with the given parameters
func GetScaleTransformMatrix3D(xscale, yscale, zscale float64) mat64.Matrix {
	return mat64.NewDense(4, 4,
		[]float64{
			xscale, 0, 0, 0,
			0, yscale, 0, 0,
			0, 0, zscale, 0,
			0, 0, 0, 1})
}

// GetScaleTransformMatrix returns a transformation matrix for a scale
// transform with the given parameters
func GetScaleTransformMatrix(xscale, yscale float64) mat64.Matrix {
	// TODO update to 3D
	return mat64.NewDense(3, 3,
		[]float64{
			xscale, 0, 0,
			0, yscale, 0,
			0, 0, 1})
}

// GetTranslationTransformMatrix3D returns a transformation matrix for a
// translation transform with the given parameters
func GetTranslationTransformMatrix3D(dx, dy, dz float64) mat64.Matrix {
	return mat64.NewDense(4, 4,
		[]float64{
			1, 0, 0, dx,
			0, 1, 0, dy,
			0, 0, 1, dz,
			0, 0, 0, 1})
}

// GetTranslationTransformMatrix returns a transformation matrix for a
// translation transform with the given parameters
func GetTranslationTransformMatrix(dx, dy float64) mat64.Matrix {
	return mat64.NewDense(3, 3,
		[]float64{
			1, 0, dx,
			0, 1, dy,
			0, 0, 1})
}

// GetRotationTransformMatrix returns a transformation matrix for a
// rotation transform with the given parameters
func GetRotationTransformMatrix(deg int) mat64.Matrix {
	// TODO update to 3D
	fdeg := float64(deg)
	frad := fdeg * (math.Pi / 180)
	return mat64.NewDense(3, 3,
		[]float64{
			math.Cos(frad), -math.Sin(frad), 0,
			math.Sin(frad), math.Cos(frad), 0,
			0, 0, 1})
}

// GetPerspectiveProjectionMatrix returns the matrix used to transform 3D
// vertices to the canonical view volume for a perspective projection
func GetPerspectiveProjectionMatrix(PRP *mat64.Vector, ViewPlane Port2D,
	back float64) mat64.Matrix {

	polygons := mat64.NewDense(4, 4, nil)
	PRPx := PRP.At(0, 0)
	PRPy := PRP.At(1, 0)
	PRPz := PRP.At(2, 0)
	umin := ViewPlane.XMin
	umax := ViewPlane.XMax
	vmin := ViewPlane.YMin
	vmax := ViewPlane.YMax

	// Row 1
	polygons.Set(0, 0,
		((2 * PRPz) / ((umax - umin) * (PRPz - back))))

	polygons.Set(0, 1, 0)

	polygons.Set(0, 2,
		(((umax + umin) - (2 * PRPx)) / ((umax - umin) * (PRPz - back))))

	polygons.Set(0, 3,
		-(((umax + umin) * PRPz) / ((umax - umin) * (PRPz - back))))

	// Row 2
	polygons.Set(1, 0, 0)

	polygons.Set(1, 1,
		((2 * PRPz) / ((vmax - vmin) * (PRPz - back))))

	polygons.Set(1, 2,
		(((vmax + vmin) - (2 * PRPy)) / ((vmax - vmin) * (PRPz - back))))

	polygons.Set(1, 3,
		-(((vmax + vmin) * PRPz) / ((vmax - vmin) * (PRPz - back))))

	// Row 3
	polygons.Set(2, 0, 0)

	polygons.Set(2, 1, 0)

	polygons.Set(2, 2,
		(1 / (PRPz - back)))

	polygons.Set(2, 3,
		-(PRPz / (PRPz - back)))

	// Row 4
	polygons.SetRow(3, []float64{0, 0, 0, 1})

	return polygons
}

// GetOrthographicProjectionMatrix returns the matrix used to transform 3D
// vertices to the canonical view volume for a orthographic projection
func GetOrthographicProjectionMatrix(PRP *mat64.Vector, ViewPlane Port2D, front,
	back float64) mat64.Matrix {

	polygons := mat64.NewDense(4, 4, nil)
	PRPx := PRP.At(0, 0)
	PRPy := PRP.At(1, 0)
	PRPz := PRP.At(2, 0)
	umin := ViewPlane.XMin
	umax := ViewPlane.XMax
	vmin := ViewPlane.YMin
	vmax := ViewPlane.YMax

	// Row 1
	polygons.Set(0, 0,
		(2 / (umax - umin)))

	polygons.Set(0, 1, 0)

	polygons.Set(0, 2,
		(((umax + umin) - (2 * PRPx)) / ((umax - umin) * PRPz)))

	polygons.Set(0, 3,
		-((umax + umin) / 2))

	// Row 2
	polygons.Set(1, 0, 0)

	polygons.Set(1, 1,
		(2 / (vmax - vmin)))

	polygons.Set(1, 2,
		(((vmax + vmin) - (2 * PRPy)) / ((vmax - vmin) * PRPz)))

	polygons.Set(1, 3,
		-((vmax + vmin) / 2))

	// Row 3
	polygons.Set(2, 0, 0)

	polygons.Set(2, 1, 0)

	polygons.Set(2, 2,
		(1 / (front - back)))

	polygons.Set(2, 3,
		-(front / (front - back)))

	// Row 4
	polygons.SetRow(3, []float64{0, 0, 0, 1})

	return polygons
}

// Normalize3Vector normalizes a 3D vector
func Normalize3Vector(vector *mat64.Vector) {
	normFactor := math.Sqrt(math.Pow(vector.At(0, 0), 2) +
		math.Pow(vector.At(1, 0), 2) +
		math.Pow(vector.At(2, 0), 2))

	vector.SetVec(0, vector.At(0, 0)/normFactor)
	vector.SetVec(1, vector.At(1, 0)/normFactor)
	vector.SetVec(2, vector.At(2, 0)/normFactor)

}

// GetViewMatrix returns the matrix to be used to transform vertices to a camera
// view based on the given parameters.
func GetViewMatrix(VPN, VUP, VRP *mat64.Vector) mat64.Matrix {
	// Calculate the camera coordinate axes
	n := VPN
	u := VectorCrossProduct(VUP, VPN)
	v := VectorCrossProduct(n, u)
	Normalize3Vector(n)
	Normalize3Vector(u)
	Normalize3Vector(v)

	// Get a translation matrix of -VRP
	translationMatrix :=
		GetTranslationTransformMatrix3D(-VRP.At(0, 0), -VRP.At(1, 0), -VRP.At(2, 0))

	// Create a rotation matrix based on the camera axes
	rotationMatrix := mat64.NewDense(4, 4, []float64{
		u.At(0, 0), u.At(1, 0), u.At(2, 0), 0,
		v.At(0, 0), v.At(1, 0), v.At(2, 0), 0,
		n.At(0, 0), n.At(1, 0), n.At(2, 0), 0,
		0, 0, 0, 1})
	viewMatrix := mat64.NewDense(4, 4, nil)
	viewMatrix.Mul(rotationMatrix, translationMatrix)

	return viewMatrix
}

// VectorCrossProduct returns the cross product of two 3d vectors
func VectorCrossProduct(a, b *mat64.Vector) *mat64.Vector {
	result := mat64.NewVector(3, []float64{

		a.At(1, 0)*b.At(2, 0) - a.At(2, 0)*b.At(1, 0),
		a.At(2, 0)*b.At(0, 0) - a.At(0, 0)*b.At(2, 0),
		a.At(0, 0)*b.At(1, 0) - a.At(1, 0)*b.At(0, 0)})

	return result

}

func setCol(mat *mat64.Dense, j int, vec *mat64.Vector) {
	mat.SetCol(j, stripVector(vec))
}

// StripVector returns a slice of the values in a mat64.Vector
func stripVector(vector *mat64.Vector) []float64 {
	var res []float64
	for i := 0; i < vector.Len(); i++ {
		res = append(res, vector.At(i, 0))
	}

	return res
}

/* Line */

// Line is a line
type Line struct {
	a, b *Vertex
}

// CreateLine creates a line with the specified vertices and
// a default clipping algorithm
func CreateLine(a, b *Vertex) *Line {
	//TODO use this where we make lines
	polygons := &Line{}
	polygons.a = a
	polygons.b = b
	return polygons
}

// Draw draws the line to the buffer
func (line *Line) Draw(buffer *SoftFrameBuffer) {

	if line.handleEdgeCases(buffer) {
		return
	}

	// Unpack the line
	qx, qy, rx, ry := line.UnpackInt()

	var x, y, start, end int
	var dy, dx, fx, fy float64

	white := RGBColor{255, 255, 255}

	m := float64(ry-qy) / float64(rx-qx)

	if math.Abs(m) < 1 {
		dx = 1.0
		dy = m
		if qx <= rx {
			x = qx
			y = qy
			start = qx
			end = rx
		} else {
			x = rx
			y = ry
			start = rx
			end = qx
		}

	} else {
		dx = 1 / m
		dy = 1.0
		if qy <= ry {
			x = qx
			y = qy
			start = qy
			end = ry
		} else {
			x = rx
			y = ry
			start = ry
			end = qy
		}
	}

	fx = float64(x)
	fy = float64(y)
	for i := start; i <= end; i++ {
		buffer.WritePixel(round(fx), round(fy), white)

		fx = fx + dx
		fy = fy + dy
	}

}

// Transform applies the transformation matrix to the line
func (line *Line) Transform(matrix mat64.Matrix) {
	// TODO can't transform lines right now
}

// Discretize rounds all of the points in the line to integers
func (line *Line) Discretize() {
	// TODO lines can't be discretized right now
}

// Project projects the line from 3D to 2D
func (line *Line) Project(projectionMatrix mat64.Matrix, orthographic bool, scene *Scene) bool {
	//TODO lines can't be projected right now
	return false
}

// Print prints the points in this line for debugging
func (line *Line) String() string {
	x0, y0, z0, x1, y1, z1 := line.UnpackFloat()
	return fmt.Sprintf("[%f,%f,%f]\n[%f,%f,%f]\n", x0, y0, z0, x1, y1, z1)
}

func round(f float64) int {
	dec := f - math.Floor(f)
	if dec >= 0.5 {
		return int(f + 1)
	}
	return int(f)
}

func findYC(x0, x1, y0, y1, clip float64) float64 {
	return ((clip-x0)/(x1-x0))*(y1-y0) + y0
}

func findXC(x0, x1, y0, y1, clip float64) float64 {
	return ((y1-clip)/(y1-y0))*(x0-x1) + x1
}

// UnpackInt returns this line's vertex coordinates as a sequence of integers
func (line *Line) UnpackInt() (qx, qy, rx, ry int) {
	qx = round(line.a.attributes[0])
	qy = round(line.a.attributes[1])
	rx = round(line.b.attributes[0])
	ry = round(line.b.attributes[1])

	return
}

// UnpackFloat returns this line's vertex coordinates as a sequence of floats
func (line *Line) UnpackFloat() (qx, qy, qz, rx, ry, rz float64) {
	qx = line.a.attributes[0]
	qy = line.a.attributes[1]
	qz = line.a.attributes[2]
	rx = line.b.attributes[0]
	ry = line.b.attributes[1]
	rz = line.b.attributes[2]

	return
}

func (line *Line) handleEdgeCases(buffer *SoftFrameBuffer) bool {
	var length int
	var variable *int
	handled := false
	black := RGBColor{0, 0, 0}

	qx, qy, rx, ry := line.UnpackInt()

	// Just a single point
	if qx == rx && qy == ry {
		buffer.WritePixel(qx, qy, black)
		return true
	}

	if qy == ry { // Horizontal
		// Check if we need to swap the points
		if rx < qx {
			swap(&rx, &qx)
			swap(&ry, &qy)
		}
		variable = &qx
		length = rx - qx
		handled = true
	} else if qx == rx { // Vertical
		if ry < qy {
			swap(&rx, &qx)
			swap(&ry, &qy)
		}
		variable = &qy
		length = ry - qy
		handled = true
	}

	for i := 0; i < length; i++ {
		buffer.WritePixel(qx, qy, black)
		(*variable)++
	}

	return handled

}

func swap(a, b *int) {
	temp := *a
	*a = *b
	*b = temp
}

/* Geometry */

// Geometry is some geometric form, comprised of vertices
type Geometry struct {
	vertices  []*Vertex
	baseColor RGBColor
}

// Draw draws this Geometry to the provided buffer
func (geo Geometry) Draw(buffer *SoftFrameBuffer) {

	// Draw the lines formws by the points of the geometry

	for i := 0; i < len(geo.vertices)-1; i++ {
		line := &Line{
			a: geo.vertices[i],
			b: geo.vertices[i+1]}

		//fmt.Printf("%s\n", geo.vertices[i].String())
		line.Draw(buffer)
	}

	// TODO adapt for 3D
	// Fill the geometry with black
	//geo.scanFill(buffer, RGBColor{0, 0, 0})

}

// Clip clips the geometry to the given port
func (geo Geometry) Clip(port Port2D) bool {
	//geo.clippingFunc(geo, port)
	return true
}

// Transform applies the given transformation matrix to the
// geometry
func (vertex *Vertex) Transform(matrix mat64.Matrix) {
	// Create a column vector from the vertex
	vector := mat64.NewVector(len(vertex.attributes), vertex.attributes)
	//fmt.Printf("%v\n", mat64.Formatted(vector))

	// Apply the transformation
	vector.MulVec(matrix, vector)

	// Put the new vertex data back
	vertex.attributes = stripVector(vector)
}

func getPerspectiveBitCode(vertex *Vertex, scene *Scene) (code int) {
	PRPz := scene.PRP.At(2, 0)
	front := scene.FrontClippingPlane
	back := scene.BackClippingPlane

	zmin := (PRPz - front) / (back - PRPz)

	x := vertex.attributes[0]
	y := vertex.attributes[1]
	z := vertex.attributes[2]

	aboveC := 1
	belowC := 1 << 1
	rightC := 1 << 2
	leftC := 1 << 3
	behindC := 1 << 4
	frontC := 1 << 5

	if y > -z {
		code |= aboveC
	}

	if y < z {
		code |= belowC
	}

	if x > -z {
		code |= rightC
	}

	if x < z {
		code |= leftC
	}

	if z < -1 {
		code |= behindC
	}

	if z > zmin {
		code |= frontC
	}

	return

}

func getOrthogonalBitCode(vertex *Vertex, scene *Scene) (code int) {

	x := vertex.attributes[0]
	y := vertex.attributes[1]
	z := vertex.attributes[2]

	aboveC := 1
	belowC := 1 << 1
	rightC := 1 << 2
	leftC := 1 << 3
	behindC := 1 << 4
	frontC := 1 << 5

	if y > 1 {
		code |= aboveC
	}

	if y < -1 {
		code |= belowC
	}

	if x > 1 {
		code |= rightC
	}

	if x < -1 {
		code |= leftC
	}

	if z < -1 {
		code |= behindC
	}

	if z > 0 {
		code |= frontC
	}

	return

}

// Discretize rounds all vertices of the geometry to integers
func (geo Geometry) Discretize() {
	for i, vertex := range geo.vertices {
		// only round x and y
		for j := 0; j < 2; j++ {
			rounded := round(vertex.attributes[j])
			geo.vertices[i].attributes[j] = float64(rounded)
		}
	}
}

func (geo Geometry) sutherlandHodgemanClip(port Port2D) {

	// Our input points for this clipping edge
	v := geo.vertices
	vprime := []*Vertex{}

	tl := Create2DVertex(port.XMin, port.YMax)

	tr := Create2DVertex(port.XMax, port.YMax)

	bl := Create2DVertex(port.XMin, port.YMin)

	br := Create2DVertex(port.XMax, port.YMin)

	corners := []*Vertex{tl, tr, br, bl}

	northInside := func(vertex *Vertex) bool {
		return vertex.attributes[1] <= float64(port.YMax)
	}

	eastInside := func(vertex *Vertex) bool {
		return vertex.attributes[0] <= float64(port.XMax)
	}

	southInside := func(vertex *Vertex) bool {
		return vertex.attributes[1] >= float64(port.YMin)
	}

	westInside := func(vertex *Vertex) bool {
		return vertex.attributes[0] >= float64(port.XMin)
	}

	insideCheckers := []func(*Vertex) bool{northInside, eastInside,
		southInside, westInside}

	// Go through each clipping edge
	cornerIdx := 0
	var clippingEdge *Line
	for _, checker := range insideCheckers {
		// Set the clipping edge line
		if cornerIdx == 3 { // Last side (west)
			clippingEdge = &Line{a: corners[cornerIdx], b: corners[0]} // Wrap
		} else {
			clippingEdge = &Line{a: corners[cornerIdx], b: corners[cornerIdx+1]}
		}

		// Check if each vertex is inside the edge and clip it to its side's
		// intersection if not
		//fmt.Printf("%v\n", v)
		var lastAppendedVertex *Vertex
		for i := 0; i < len(v)-1; i++ {
			side := &Line{a: v[i], b: v[i+1]}
			intersection := findIntersection(side, clippingEdge)

			if checker(side.a) && !checker(side.b) {
				if !side.a.Equal(lastAppendedVertex) {
					vprime = append(vprime, side.a)
				}
				vprime = append(vprime, intersection)
				lastAppendedVertex = intersection
			} else if !checker(side.a) && checker(side.b) {
				vprime = append(vprime, intersection)
				vprime = append(vprime, side.b)
				lastAppendedVertex = side.b
			} else if checker(side.a) && checker(side.b) {
				if !side.a.Equal(lastAppendedVertex) {
					vprime = append(vprime, side.a)
				}
				vprime = append(vprime, side.b)
				lastAppendedVertex = side.b
			}
		}

		// Remove duplicate vertices
		removeDuplicateVertexPtrs(vprime)

		// Make sure this clipped polygon is closed
		if len(vprime) >= 1 && !vprime[0].Equal(vprime[len(vprime)-1]) {
			vprime = append(vprime,
				Create2DVertex(vprime[0].attributes[0],
					vprime[0].attributes[1]))

		}

		cornerIdx++
		v = make([]*Vertex, len(vprime))
		copy(v, vprime)
		vprime = make([]*Vertex, 0)
	}

	// Update the polygon
	geo.vertices = v

}

func findIntersection(source, target *Line) *Vertex {
	//var x, y int

	// First we will calculate/store a few points/vectors we need for
	// computation

	P0 := source.a
	P1 := source.b
	T0 := target.a
	T1 := target.b

	// We need the vector form of our source line (P1-P0)
	D := VAdd(P1, VScale(P0, -1))
	T := VAdd(T1, VScale(T0, -1))

	// We need an endpoint of our target line
	Pe := target.a

	// And we need the normal vector to our target line
	N := T.TwoDNormal()

	// We use the parametric form of our source line, P(t) = P0 + tD
	// and assume P(t) is on the target line.  This would have to satisfy
	// N dot [P(t) - Pe] = 0, since the dot product of a segment and its
	// normal is always zero.  We can solve for the parameter t for the
	// intersection, if any.  The line actually intersects iff 0 <= t <= 1,
	// in which case we plug t back into P(t) to get the intersection point.
	// Yay math!  The formula we use to get t is just an algebraic manipulation
	// of the above equation to solve for t.

	M := VAdd(P0, VScale(Pe, -1))

	t := -1 * (VDot(N, M) / VDot(N, D))

	// P(t) = P0 + tD

	if t >= 0 && t <= 1 {
		return VAdd(P0, VScale(D, t))
	}
	return nil

}

// AddVertex adds a vertex to the geometry
func (geo Geometry) AddVertex(vertex *Vertex) {
	geo.vertices = append(geo.vertices, vertex)
}

// Print prints the points in this Geometry
func (geo Geometry) Print() {
	for _, vertex := range geo.vertices {
		fmt.Printf("%s\n", vertex.String())
	}
}

// scanFill fills the Geometry with the scan technique.  Only meant for 2D
// geometry
func (geo Geometry) scanFill(buffer *SoftFrameBuffer, scene *Scene) {
	extremas := &XSortedVertexSlice{}
	// Scan every row in the buffer
	for lineY := 0; lineY < buffer.Height; lineY++ {
		scanLine := CreateLine(Create3DVertex(0, float64(lineY), 0),
			Create3DVertex(float64(buffer.Width), float64(lineY), 0))

		// Clear the extrema list
		extremas.values = nil

		// Go through each edge
		for vertIdx := 0; vertIdx < len(geo.vertices)-1; vertIdx++ {
			edge := CreateLine(geo.vertices[vertIdx], geo.vertices[vertIdx+1])
			if !isEdgeValid(edge, lineY) {
				continue
			}

			// Find the intersection point of the scan line and this edge
			intersection := findIntersection(edge, scanLine)
			/*
				if edge.a.attributes[2] != 0 ||
					edge.b.attributes[2] != 0 {
					fmt.Fprintf(os.Stderr, "%s\n", edge.String())
				}
			*/

			// No intersection
			if intersection == nil {
				continue
			}
			// Interpolate the z value for the intersection
			interpolateEdgeZ(edge.a, edge.b, intersection)

			extremas.values =
				append(extremas.values, intersection)
		}

		// No extremas for this scan line
		if len(extremas.values) == 0 {
			continue
		}

		// Sort extremas by X value
		sort.Sort(extremas)

		// Scan over the line
		fill := false
		extremeIdx := 0

		for x := 0; x < buffer.Width; x++ {
			if fill {
				interpolatedPt := Create3DVertex(float64(x), float64(lineY), 0)
				// Interpolate a Z value for this point if needed
				if extremeIdx >= 1 {
					interpolateScanlineZ(extremas.values[extremeIdx], extremas.values[extremeIdx-1],
						interpolatedPt)
				}
				z := interpolatedPt.attributes[2]
				//fmt.Printf("%f\n", z)
				//fmt.Fprintf(os.Stderr, "%f\n", z)

				var front, back float64
				if scene.UseOrthographicProjection {
					front = 0
					back = -1
				} else {
					front = (scene.PRP.At(2, 0) - scene.FrontClippingPlane) /
						(scene.BackClippingPlane - scene.PRP.At(2, 0))
					back = -1
				}
				// Check and update the Z buffer
				if buffer.updateZBuffer(x, lineY, z) {
					buffer.writeDepthCueColor(z, front, back, x, lineY, geo.baseColor)
				}
			}
			for extremeIdx < len(extremas.values) &&
				x == round(extremas.values[extremeIdx].attributes[0]) {
				fill = !fill
				extremeIdx++
			}
		}
	}
}

func interpolateEdgeZ(a, b, target *Vertex) {
	z1 := a.attributes[2]
	z2 := b.attributes[2]
	y1 := a.attributes[1]
	y2 := b.attributes[1]
	yTarget := target.attributes[1]

	var interpZ float64
	// Check if z value for target vertex is already set
	if target.attributes[2] != 0 {
		interpZ = target.attributes[2]
	} else {
		interpZ = z1 - (z1-z2)*((y1-yTarget)/(y1-y2))
	}
	target.attributes[2] = interpZ

}
func interpolateScanlineZ(a, b, target *Vertex) {
	za := a.attributes[2]
	zb := b.attributes[2]
	xa := a.attributes[0]
	xb := b.attributes[0]
	xp := target.attributes[0]

	var interpZ float64
	// Check if z value for target vertex is already set
	if target.attributes[2] != 0 {
		interpZ = target.attributes[2]
	} else {
		interpZ = zb - (zb-za)*((xb-xp)/(xb-xa))
	}
	target.attributes[2] = interpZ

}

func isEdgeValid(edge *Line, scanY int) bool {
	// Is the edge horizontal?
	if edge.a.attributes[1] == edge.b.attributes[1] {
		return false
	}

	// Is the top vertex of the edge on the scan line?
	if int(maxY(edge.a, edge.b).attributes[1]) == scanY {
		return false
	}

	return true
}

/* XSortedVertex */

// XSortedVertexSlice is a wrapper for when we want to sort Vertex
// slices by their x values.  Implements golangs sorting interface
type XSortedVertexSlice struct {
	values []*Vertex
}

// Len returns the length of the slice.  For the sort interface
func (slice *XSortedVertexSlice) Len() int {
	return len(slice.values)
}

// Less returns whether Vertex at index i is less than that at index j.
// For the sort interface
func (slice *XSortedVertexSlice) Less(i, j int) bool {
	return slice.values[i].attributes[0] < slice.values[j].attributes[0]
}

// Swap swaps the Vertex at index i with the one at index j.
// For the sort interface
func (slice *XSortedVertexSlice) Swap(i, j int) {
	temp := slice.values[i]
	slice.values[i] = slice.values[j]
	slice.values[j] = temp
}

/*
func parsePolygonObject(lines []string) ([]*Geometry, error) {
	var res []*Geometry
	toAdd := &Geometry{}
	toAdd.clippingFunc = (*Geometry).sutherlandHodgemanClip
	polygonFinished := true
	failed := false
	for _, line := range lines {
		line = strings.Trim(line, " ")
		tokens := strings.Fields(line)
		delim := tokens[len(tokens)-1]

		// A new polygon must begin with a moveto command
		if polygonFinished && delim != "moveto" {
			failed = true
			break
		}

		if delim == "stroke" {
			// Final point and first point must be the same
			if !toAdd.vertices[0].Equal(toAdd.vertices[len(toAdd.vertices)-1]) {
				return nil, errors.New("Detected unclosed polygon")
			}
			res = append(res, toAdd)
			polygonFinished = true
			toAdd = &Geometry{}
			toAdd.vertices = nil
			toAdd.clippingFunc = (*Geometry).sutherlandHodgemanClip
			continue
		}

		polygonFinished = false

		if !(delim == "moveto" || delim == "lineto") {
			failed = true
			break
		}

		coordStrs := tokens[:len(tokens)-1]

		// We need exactly two coordinates in every
		// polygon line
		if len(coordStrs) != 2 {
			failed = true
			break
		}

		// Add the current point
		xCoord, err := strconv.Atoi(coordStrs[0])
		if err != nil {
			failed = true
		}

		yCoord, err := strconv.Atoi(coordStrs[1])
		if err != nil {
			failed = true
		}

		point := Create2DVertex(float64(xCoord), float64(yCoord))

		toAdd.AddVertex(point)

	}

	if failed {
		return nil, errors.New("Failed parsing a polygon\n")
	}

	if !polygonFinished {
		return nil, errors.New("Detected incomplete polygon specification\n")
	}

	return res, nil
}
*/

// File represents an input file that can be parsed for different input objects
type File struct {
	filePath string
	handle   *os.File
	scanner  *bufio.Scanner
	lineIdx  int
}

// SMFFile is an SMF input file
type SMFFile struct {
	Fp                     *File
	vertexDelim, faceDelim string
	BaseColor              RGBColor
}

// OpenSMFFile opens the SMF file at the given path
func OpenSMFFile(path string) (*SMFFile, error) {

	fp, err := os.Open(path)
	if err != nil {
		return nil, err
	}

	scanner := bufio.NewScanner(fp)

	file := &File{
		filePath: path,
		handle:   fp,
		scanner:  scanner,
		lineIdx:  0}

	smfFile := &SMFFile{
		Fp: file}

	// Supported line tokens
	smfFile.vertexDelim = "v"
	smfFile.faceDelim = "f"

	return smfFile, err
}

// Close closes this SMFFile
func (file *File) Close() {
	file.handle.Close()
}

// ParseSMFObjects parses objects from an SMF file
func (file *SMFFile) ParseSMFObjects() ([]Geometry, []*Vertex, error) {
	scanner := file.Fp.scanner
	var err error
	vertices := []*Vertex{}
	var polygons []Geometry
	for scanner.Scan() {
		line := scanner.Text()

		// Ignore blank lines
		if strings.Trim(line, " ") == "" {
			continue
		}
		tokens := strings.Fields(line)

		if len(tokens) < 1 {
			err = fmt.Errorf("Bad SMF line format on line %d",
				file.Fp.lineIdx)
			break
		}
		if tokens[0] == file.vertexDelim {
			if len(tokens) < 4 {
				err = fmt.Errorf("Line %d: Three coordinates required", file.Fp.lineIdx)
				break
			}
			var x, y, z float64
			x, err = strconv.ParseFloat(tokens[1], 64)
			y, err = strconv.ParseFloat(tokens[2], 64)
			z, err = strconv.ParseFloat(tokens[3], 64)

			if err != nil {
				break
			}

			vertex := Create3DVertex(x, y, z)
			vertices = append(vertices, vertex)

		} else if tokens[0] == file.faceDelim {
			geo := Geometry{}
			if len(tokens) < 4 {
				err = fmt.Errorf("Line %d: A face needs at least three"+
					" vertices", file.Fp.lineIdx)
				break
			}
			for i := 1; i < len(tokens); i++ {
				var vidx int
				vidx, err = strconv.Atoi(tokens[i])
				if err != nil {
					break
				}
				if vidx < 1 || vidx > len(vertices) {
					err = fmt.Errorf("Line %d: Invalid vertex index", file.Fp.lineIdx)
					break
				}
				geo.vertices = append(geo.vertices, vertices[vidx-1]) // index starts at 1 in file
			}

			if err != nil {
				break
			}

			if len(geo.vertices) < 1 {
				fmt.Println()
			}

			// Add the "closing" vertex
			closing := geo.vertices[0]
			geo.vertices = append(geo.vertices, closing)

			// Set the base color
			geo.baseColor = file.BaseColor

			polygons = append(polygons, geo)
			//geo.Print()

		} else {
			err = fmt.Errorf("Line %d: Invalid SMF token", file.Fp.lineIdx)
			break
		}

		file.Fp.lineIdx++
	}

	return polygons, vertices, err
}

/* XPMFile */

// XPMFile is an XPM image file
type XPMFile struct {
	XPMFileFormatStr string

	filePath string
	handle   *os.File
	writer   *bufio.Writer
}

// OpenXPMFile opens a new XPM file for writing
func OpenXPMFile(path string) (*XPMFile, error) {
	//fp, err := os.Create(path)
	fp := os.Stdout
	err := errors.New("")

	writer := bufio.NewWriter(fp)

	file := &XPMFile{
		filePath: path,
		handle:   fp,
		writer:   writer}

	file.XPMFileFormatStr =
		"/* XPM */\n" +
			"static char *sco100[] = {\n" +
			"/* width height num_colors chars_per_pixel */\n" +
			"\"%d %d %d %d\",\n" +
			"/* colors */\n" +
			"%s" +
			"/* pixels */\n" +
			"%s" +
			"};\n"

	return file, err
}

// WriteSoftBuffer writes the given buffer to an xpm file
func (file *XPMFile) WriteSoftBuffer(buffer *SoftFrameBuffer) error {
	// Build the color list string
	var colorList string
	for k, v := range buffer.colors {
		colorList += fmt.Sprintf("\"%s c #%s\",\n", v, k.GenerateHexString())
	}

	// Build the rows string
	var rowsList bytes.Buffer
	for y := 0; y < buffer.Height; y++ {
		rowsList.WriteString("\"")
		for x := 0; x < buffer.Width; x++ {
			rowsList.WriteString(buffer.colors[buffer.buffer[y][x]])
		}
		rowsList.WriteString("\"")
		if y != buffer.Height-1 {
			rowsList.WriteString(",")
		}
		rowsList.WriteString("\n")
	}

	// Put it all together
	fileStr := fmt.Sprintf(file.XPMFileFormatStr, buffer.Width,
		buffer.Height, len(buffer.colors), 3, colorList, rowsList.String())

	_, err := file.writer.WriteString(fileStr)
	file.writer.Flush()

	return err
}

/* ----Misc routines---- */

func getCwd() string {
	dir, err := filepath.Abs(filepath.Dir(os.Args[0]))
	if err != nil {
		log.Fatal(err)
	}
	return dir
}

// GenerateHexString returns the hex string for this color
func (color RGBColor) GenerateHexString() string {
	return fmt.Sprintf("%02x%02x%02x", color.R, color.G, color.B)
}

func minX(a, b *Vertex) *Vertex {
	if a.attributes[0] < b.attributes[0] {
		return a
	}

	return b
}

func minY(a, b *Vertex) *Vertex {
	if a.attributes[1] < b.attributes[1] {
		return a
	}

	return b
}

func maxX(a, b *Vertex) *Vertex {
	if a.attributes[0] >= b.attributes[0] {
		return a
	}

	return b
}

func maxY(a, b *Vertex) *Vertex {
	if a.attributes[1] >= b.attributes[1] {
		return a
	}

	return b
}

/*
func removeDuplicateVertexPtrs(slice []*Vertex) []*Vertex {
	found := make(map[*Vertex]bool)
	res := make([]*Vertex, len(slice))
	originalLength := len(slice)
	newLength := 0

	for i := 0; i < originalLength; i++ {
		vertex := slice[i]
		_,included := found[vertex]
		if !included {
			res[newLength] = vertex
			newLength++
		}
		found[vertex] = true
	}

	return res[:newLength]
}
*/

func removeDuplicateVertexPtrs(slice []*Vertex) {
	found := make(map[*Vertex]bool)
	length := len(slice)

	for i := 0; i < length; i++ {
		vertex := slice[i]
		_, included := found[vertex]
		if included {
			slice = append(slice[:i], slice[i+1:]...)
			length--
			i--
		}
		found[vertex] = true
	}
}
