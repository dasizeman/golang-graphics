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
)

/* ----Interfaces---- */

// Drawable is some kind of geometry that can be rasterized to a frame buffer
type Drawable interface {
	Draw(buffer *SoftFrameBuffer)
	clip(buffer *SoftFrameBuffer) bool
}

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
	colors        map[RGBColor]string
}

// NewSoftFrameBuffer allocates a new SoftFrameBuffer with specified width and
// height, with 8bit color channels
func NewSoftFrameBuffer(width, height int) *SoftFrameBuffer {
	newBuffer := &SoftFrameBuffer{
		Width:  width,
		Height: height}

	newBuffer.colors = make(map[RGBColor]string)

	white := RGBColor{255, 255, 255}
	black := RGBColor{0, 0, 0}
	// Hard code black and white for now
	newBuffer.colors[white] = "#"
	newBuffer.colors[black] = "*"

	// Allocate the pixel buffer
	newBuffer.buffer = make([][]RGBColor, newBuffer.Height+1)
	for y := range newBuffer.buffer {
		newBuffer.buffer[y] = make([]RGBColor, newBuffer.Width+1)
	}

	// Set the buffer to white
	for y := 0; y <= newBuffer.Height; y++ {
		for x := 0; x <= newBuffer.Width; x++ {
			newBuffer.WritePixel(x, y, white)
		}
	}

	return newBuffer
}

// WritePixel writes a pixel color to the buffer, assuming (0,0) is the
//bottom-left corner
func (frameBuffer *SoftFrameBuffer) WritePixel(x, y int, color RGBColor) {

	// We will reverse the y origin for compatibility with XPM output
	frameBuffer.buffer[frameBuffer.Height-y][x] = color

	// Add this to our set of colors
	//frameBuffer.colors[color] = GenerateHexString(color)

}

/* Scene */

// Scene is a virtual space that can be rendered to a frame buffer
type Scene struct {
	Objects []Drawable
}

// Render renders the scene to the given buffer
func (scene *Scene) Render(buffer *SoftFrameBuffer) {
	for _, object := range scene.Objects {
		object.Draw(buffer)
	}
}

/* Vertex */

// Vertex is a vertex with attributes
type Vertex struct {
	attributes []float64
}

// Create2DVertex creates a vertex with the specified x and y coordinates
func Create2DVertex(a, b float64) *Vertex {
	// TODO use this where we create 2D vertices
	result := &Vertex{}
	result.AddAttribute(a)
	result.AddAttribute(b)
	return result
}

// AddAttribute adds an attribute to vertex
func (vertex *Vertex) AddAttribute(attr float64) {
	vertex.attributes = append(vertex.attributes, attr)
}

func (vertex *Vertex) String() string {
	return fmt.Sprintf("%v", vertex.attributes)
}

// Equal checks if this Vertex is equal to another Vertex
func (vertex *Vertex) Equal(other *Vertex) bool {
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

// Print prints this vertex
func (vertex *Vertex) Print() {
	fmt.Printf("(")
	for _, attr := range vertex.attributes {
		fmt.Printf("%f,", attr)
	}
	fmt.Printf(")\n")
}

/* Line */

// Line is a line
type Line struct {
	a, b              *Vertex
	clippingAlgorithm func(*Line, int, int) bool
}

// CreateLine creates a line with the specified vertices and
// a default clipping algorithm
func CreateLine(a, b *Vertex) *Line {
	//TODO use this where we make lines
	result := &Line{}
	result.a = a
	result.b = b
	result.clippingAlgorithm = (*Line).cohenSutherlandClip
	return result
}

func (line *Line) clip(buffer *SoftFrameBuffer) bool {
	return line.clippingAlgorithm(line, buffer.Width, buffer.Height)
}

// Print prints the points in this line for debugging
func (line *Line) Print() {
	x0, y0, x1, y1 := line.UnpackFloat()
	fmt.Printf("[%f,%f]\n[%f,%f]\n", x0, y0, x1, y1)
}

// Draw draws the line to the buffer
func (line *Line) Draw(buffer *SoftFrameBuffer) {

	// Our clipper has reported that the line is completely outside the buffer
	if !line.clip(buffer) {
		return
	}

	if line.handleEdgeCases(buffer) {
		return
	}

	// Unpack the line
	qx, qy, rx, ry := line.UnpackInt()

	var x, y, start, end int
	var dy, dx, fx, fy float64

	black := RGBColor{0, 0, 0}

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
		buffer.WritePixel(round(fx), round(fy), black)

		fx = fx + dx
		fy = fy + dy
	}

}

func round(f float64) int {
	dec := f - math.Floor(f)
	if dec >= 0.5 {
		return int(f + 1)
	}
	return int(f)
}

const (
	csWest  int = 1
	csEast  int = 2
	csSouth int = 4
	csNorth int = 8
)

func (line *Line) cohenSutherlandClip(width, height int) bool {
	aCode := line.a.getCSBitcode(width, height)
	bCode := line.b.getCSBitcode(width, height)
	x0, y0, x1, y1 := line.UnpackFloat()

	// Line completely out
	if aCode&bCode != 0 {
		return false
	}

	for aCode|bCode != 0 {
		// Find the first bit flip
		bitIdx := 1
		for aCode&bitIdx == bCode&bitIdx {
			bitIdx = bitIdx << 1
		}

		// "Fix" this flip by clipping a coordinate
		var clippedVertex *Vertex
		var xc, yc int
		switch bitIdx {
		case csWest:
			xc = 0
			yc = round(findYC(x0, x1, y0, y1, 0))
			clippedVertex = minX(line.a, line.b)

		case csEast:
			xc = width
			yc = round(findYC(x0, x1, y0, y1, float64(width)))
			clippedVertex = maxX(line.a, line.b)

		case csSouth:
			yc = 0
			xc = round(findXC(x0, x1, y0, y1, 0))
			clippedVertex = minY(line.a, line.b)

		case csNorth:
			yc = height
			xc = round(findXC(x0, x1, y0, y1, float64(height)))
			clippedVertex = maxY(line.a, line.b)
		}

		clippedVertex.attributes[0] = float64(xc)
		clippedVertex.attributes[1] = float64(yc)

		// Recalculate bitcodes
		aCode = line.a.getCSBitcode(width, height)
		bCode = line.b.getCSBitcode(width, height)
	}

	return true
}

func findYC(x0, x1, y0, y1, clip float64) float64 {
	return ((clip-x0)/(x1-x0))*(y1-y0) + y0
}

func findXC(x0, x1, y0, y1, clip float64) float64 {
	return ((y1-clip)/(y1-y0))*(x0-x1) + x1
}

func (vertex *Vertex) getCSBitcode(width, height int) int {

	code := 0

	x := vertex.attributes[0]
	y := vertex.attributes[1]

	if x < 0 {
		code |= csWest
	} else if x > float64(width) {
		code |= csEast
	}

	if y < 0 {
		code |= csSouth
	} else if y > float64(height) {
		code |= csNorth
	}

	return code
}

// UnpackInt returns this line's vertex coordinates as a sequence of integers
func (line *Line) UnpackInt() (qx, qy, rx, ry int) {
	qx = int(line.a.attributes[0])
	qy = int(line.a.attributes[1])
	rx = int(line.b.attributes[0])
	ry = int(line.b.attributes[1])

	return
}

// UnpackFloat returns this line's vertex coordinates as a sequence of floats
func (line *Line) UnpackFloat() (qx, qy, rx, ry float64) {
	qx = line.a.attributes[0]
	qy = line.a.attributes[1]
	rx = line.b.attributes[0]
	ry = line.b.attributes[1]

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
	vertices     []*Vertex
	lines        []*Line
	clippingFunc func(*Geometry, int, int)
}

// Draw draws this Geometry to the provided buffer
func (geo *Geometry) Draw(buffer *SoftFrameBuffer) {
	//TODO
	geo.clip(buffer)

	// Draw the lines formws by the points of the geometry

	for i := 0; i < len(geo.vertices)-1; i++ {
		line := &Line{
			a:                 geo.vertices[i],
			b:                 geo.vertices[i+1],
			clippingAlgorithm: (*Line).cohenSutherlandClip}

		line.Draw(buffer)
	}

	/*
		for _, line := range geo.lines {
			line.Draw(buffer)
		}
	*/

	// TODO adapt for 3D
	// Fill the geometry with black
	geo.scanFill(buffer, RGBColor{0, 0, 0})

}

func (geo *Geometry) clip(buffer *SoftFrameBuffer) bool {
	geo.clippingFunc(geo, buffer.Width, buffer.Height)
	return true
}

func (geo *Geometry) sutherlandHodgemanClip(width int, height int) {

	// Our input points for this clipping edge
	v := geo.vertices
	vprime := []*Vertex{}
	fheight := float64(height)
	fwidth := float64(width)

	// TL
	tl := &Vertex{}
	tl.AddAttribute(0)
	tl.AddAttribute(fheight)

	// TR
	tr := &Vertex{}
	tr.AddAttribute(fwidth)
	tr.AddAttribute(fheight)

	// BL
	bl := &Vertex{}
	bl.AddAttribute(0)
	bl.AddAttribute(0)

	// BR
	br := &Vertex{}
	br.AddAttribute(fwidth)
	br.AddAttribute(0)

	corners := []*Vertex{tl, tr, br, bl}

	northInside := func(vertex *Vertex) bool {
		return vertex.attributes[1] <= fheight
	}

	eastInside := func(vertex *Vertex) bool {
		return vertex.attributes[0] <= fwidth
	}

	southInside := func(vertex *Vertex) bool {
		return vertex.attributes[1] >= 0
	}

	westInside := func(vertex *Vertex) bool {
		return vertex.attributes[0] >= 0
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
		for i := 0; i < len(v)-1; i++ {
			side := &Line{a: v[i], b: v[i+1]}
			intersection := findIntersection(side, clippingEdge)

			if checker(side.a) && !checker(side.b) {
				vprime = append(vprime, side.a)
				vprime = append(vprime, intersection)
			} else if !checker(side.a) && checker(side.b) {
				vprime = append(vprime, intersection)
				vprime = append(vprime, side.b)
			} else if checker(side.a) && checker(side.b) {
				vprime = append(vprime, side.a)
				vprime = append(vprime, side.b)
			}
		}
		cornerIdx++
		v = vprime
		vprime = []*Vertex{}
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
func (geo *Geometry) AddVertex(vertex *Vertex) {
	geo.vertices = append(geo.vertices, vertex)
}

// Print prints the points in this Geometry
func (geo *Geometry) Print() {
	for _, vertex := range geo.vertices {
		vertex.Print()
	}
}

// scanFill fills the Geometry with the scan technique.  Only meant for 2D
// geometry
func (geo *Geometry) scanFill(buffer *SoftFrameBuffer, color RGBColor) {
	extremas := &XSortedVertexSlice{}
	// Scan every row in the buffer
	for lineY := 1; lineY <= buffer.Height; lineY++ {
		scanLine := CreateLine(Create2DVertex(0, float64(lineY)),
			Create2DVertex(float64(buffer.Width), float64(lineY)))

		// Clear the extrema list
		extremas.values = nil

		// Go through each edge
		for _, edge := range geo.lines {
			if !isEdgeValid(edge, lineY) {
				continue
			}

			// Find the intersection point of the scan line and this edge
			intersection := findIntersection(edge, scanLine)

			// No intersection
			if intersection == nil {
				continue
			}
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
		for x := 1; x <= buffer.Width; x++ {
			if fill {
				buffer.WritePixel(x, lineY, color)
			}
			for extremeIdx < len(extremas.values) &&
				x == round(extremas.values[extremeIdx].attributes[0]) {
				fill = !fill
				extremeIdx++
			}
		}
	}
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

		point := &Vertex{}
		point.AddAttribute(float64(xCoord))
		point.AddAttribute(float64(yCoord))

		toAdd.AddVertex(point)

	}

	if failed {
		return nil, errors.New("Failed parsing a polygon\n")
	}

	if !polygonFinished {
		return nil, errors.New("Detected incomplete polygon specification\n")
	}

	// It is useful to store the edges, too
	for _, geometry := range res {
		for i := 0; i < len(geometry.vertices)-1; i++ {
			var edge *Line
			edge = CreateLine(geometry.vertices[i], geometry.vertices[i+1])
			geometry.lines = append(geometry.lines, edge)

			//TODO use these lines in polygon drawing and clipping
		}
	}

	return res, nil
}

/* PostScriptFile */

// PostScriptFile is a wrapper for a PostScript format specification of a scene
type PostScriptFile struct {
	// String "constants"
	BeginDelim, EndDelim, LineDelim string
	PolygonDelims                   map[string]bool

	filePath string
	handle   *os.File
	scanner  *bufio.Scanner
	lineIdx  int
}

// OpenPostScriptFile opens the PostScript file at the given path
func OpenPostScriptFile(path string) (*PostScriptFile, error) {

	fp, err := os.Open(path)
	if err != nil {
		return nil, err
	}

	scanner := bufio.NewScanner(fp)

	file := &PostScriptFile{
		filePath: path,
		handle:   fp,
		scanner:  scanner,
		lineIdx:  0}

	// Supported line tokens
	file.BeginDelim = "%%%BEGIN"
	file.EndDelim = "%%%END"
	file.LineDelim = "Line"
	file.PolygonDelims = make(map[string]bool)
	file.PolygonDelims["moveto"] = true
	file.PolygonDelims["lineto"] = true
	file.PolygonDelims["stroke"] = true
	file.lineIdx = 1

	return file, err
}

// ParseObjects parses all Drawable geometry objects out of a PostScriptFile
func (file *PostScriptFile) ParseObjects() ([]Drawable, error) {

	var objects []Drawable
	foundBegin := false

	scanner := file.scanner
	for !foundBegin && scanner.Scan() {
		line := scanner.Text()
		file.lineIdx++
		if line == file.BeginDelim {
			foundBegin = true
		}
	}
	if !foundBegin {
		err := errors.New("PostScript file did not have the correct BEGIN" +
			" delimiter\n")
		return objects, err
	}

	file.lineIdx++

	// This contains all of the groups of lines we find that define a Polygon
	var polygonLines []string

	for scanner.Scan() {
		line := scanner.Text()
		line = strings.Trim(line, " ")

		if line == file.EndDelim {
			break
		}

		tokens := strings.Fields(line)

		// Skip blank lines
		if line == "" {
			file.lineIdx++
			continue
		}

		delim := tokens[len(tokens)-1]

		// Polygons
		_, polygonDelimFound := file.PolygonDelims[delim]
		if polygonDelimFound {
			polygonLines = append(polygonLines, line)
			file.lineIdx++
			continue
		}

		// Single line objects
		switch delim {
		case file.LineDelim:
			object := parseLineObject(tokens[:len(tokens)-1])
			if object == nil {
				errStr :=
					fmt.Sprintf("Failed to parse %s object on line %d\n",
						file.LineDelim, file.lineIdx)
				return objects, errors.New(errStr)
			}

			objects = append(objects, object)

			/*
				// Some line debug stuff
				line := object.(*Line)
				fmt.Printf("Parsed object %#v\n", line)
				fmt.Printf("%s\n%s\n", line.a.String(), line.b.String())
			*/

		default:
			errStr := fmt.Sprintf("Unrecognized object delimiter %s on line"+
				"%d\n",
				delim, file.lineIdx)
			return objects, errors.New(errStr)

		}

		file.lineIdx++

	}

	if len(polygonLines) > 0 {
		polygons, err := parsePolygonObject(polygonLines)
		if err != nil {
			return objects, err
		}
		for _, polygon := range polygons {
			objects = append(objects, polygon)
		}
	}

	return objects, nil

}

// Close closes this PostScriptFile
func (file *PostScriptFile) Close() {
	file.handle.Close()
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
		buffer.Height, len(buffer.colors), 1, colorList, rowsList.String())

	_, err := file.writer.WriteString(fileStr)
	file.writer.Flush()

	return err
}

/* ---PostScript object parsing routines--- */

func parseLineObject(tokens []string) Drawable {
	var points []int

	for i := 0; i < 4; i++ {
		num, err := strconv.Atoi(tokens[i])
		if err != nil {
			return nil
		}
		points = append(points, num)
	}

	a := &Vertex{}
	a.AddAttribute(float64(points[0]))
	a.AddAttribute(float64(points[1]))

	b := &Vertex{}
	b.AddAttribute(float64(points[2]))
	b.AddAttribute(float64(points[3]))

	res := &Line{
		a:                 a,
		b:                 b,
		clippingAlgorithm: (*Line).cohenSutherlandClip}

	return res

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
