package graphics

import (
	"bufio"
	"errors"
	"fmt"
	"log"
	"os"
	"path/filepath"
	"strconv"
	"strings"
)

/* ----Interfaces---- */

// Drawable is some kind of geometry that can be rasterized to a frame buffer
type Drawable interface {
	Draw(buffer *SoftFrameBuffer)
	clip(buffer *SoftFrameBuffer)
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
	for y := 0; y < newBuffer.Height; y++ {
		for x := 0; x < newBuffer.Width; x++ {
			newBuffer.WritePixel(x, y, white)
		}
	}

	return newBuffer
}

// WritePixel writes a pixel color to the buffer, assuming (0,0) is the
//bottom-left corner
func (frameBuffer *SoftFrameBuffer) WritePixel(x, y int, color RGBColor) {

	// We will reverse the y origin for compatibility with XPM output
	frameBuffer.buffer[frameBuffer.Height-y-1][x] = color

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

// AddAttribute adds an attribute to vertex
func (vertex *Vertex) AddAttribute(attr float64) {
	vertex.attributes = append(vertex.attributes, attr)
}

func (vertex *Vertex) String() string {
	return fmt.Sprintf("%v", vertex.attributes)
}

/* Line */

// Line is a line
type Line struct {
	a, b *Vertex
}

func (line *Line) clip(buffer *SoftFrameBuffer) {

}

// Draw draws the line to the buffer
func (line *Line) Draw(buffer *SoftFrameBuffer) {

	line.clip(buffer)

	if line.handleEdgeCases(buffer) {
		return
	}

	// Unpack the line
	qx, qy, rx, ry := line.Unpack()

	var independentCoord, dependentCoord *int
	var decider func() int

	dx := rx - qx
	dy := ry - qy

	m := float64(dy / dx)

	decisionConstantTerm := func() int {
		return 2 * (rx*dy - ry*dx)
	}

	caseOneDecider := func() int {

	}

	// Use Bresenham's algorithm
	// Case 1: 0 < m < 1
	// Case 2: m > 1
	// Case 3: m < -1
	// Case 4: -1 < m < 0
	if m > 0 && m < 1 {

	} else if m > 1 {

	} else if m < -1 {

	} else if m > -1 && m < 0 {

	}

}

func (line *Line) Unpack() (qx, qy, rx, ry int) {
	qx = int(line.a.attributes[0])
	qy = int(line.a.attributes[1])
	rx = int(line.b.attributes[0])
	ry = int(line.b.attributes[1])

	return
}

func (line *Line) handleEdgeCases(buffer *SoftFrameBuffer) bool {
	var length int
	var variable *int
	handled := false
	black := RGBColor{0, 0, 0}

	qx, qy, rx, ry := line.Unpack()

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

/* PostScriptFile */

// PostScriptFile is a wrapper for a PostScript format specification of a scene
type PostScriptFile struct {
	// String "constants"
	BeginDelim, EndDelim, LineDelim string

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

	file.lineIdx = 1

	return file, err
}

// ParseObjects parses all Drawable geometry objects out of a PostScriptFile
func (file *PostScriptFile) ParseObjects() ([]Drawable, error) {

	var objects []Drawable
	// Check if the first line is the begin delimiter
	scanner := file.scanner
	scanner.Scan()
	line := scanner.Text()
	var err error
	if line != file.BeginDelim {
		err := errors.New("PostScript file did not have the correct BEGIN" +
			" delimiter\n")
		return objects, err
	}

	file.lineIdx++

	for scanner.Scan() {
		line = scanner.Text()

		if line == file.EndDelim {
			break
		}

		tokens := strings.Split(line, " ")

		// Skip blank lines
		if len(tokens) <= 0 {
			continue
		}

		delim := tokens[len(tokens)-1]

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
			file.lineIdx++

		default:
			errStr := fmt.Sprintf("Unrecognized object delimiter %s on line"+
				"%d\n",
				delim, file.lineIdx)
			return objects, errors.New(errStr)

		}

	}

	return objects, err

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

func OpenXPMFile(path string) (*XPMFile, error) {
	fp, err := os.Create(path)
	if err != nil {
		return nil, err
	}

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
		colorList += fmt.Sprintf("\"%s c #%s\",\n", v, GenerateHexString(k))
	}

	// Build the rows string
	var rowsList string
	for y := 0; y < buffer.Height; y++ {
		rowsList += "\""
		for x := 0; x < buffer.Width; x++ {
			rowsList += buffer.colors[buffer.buffer[y][x]]
		}
		rowsList += "\""
		if y != buffer.Height-1 {
			rowsList += ","
		}
		rowsList += "\n"
	}

	// Put it all together
	fileStr := fmt.Sprintf(file.XPMFileFormatStr, buffer.Width,
		buffer.Height, len(buffer.colors), 1, colorList, rowsList)

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

	return &Line{a, b}

}

/* ----Misc routines---- */

func getCwd() string {
	dir, err := filepath.Abs(filepath.Dir(os.Args[0]))
	if err != nil {
		log.Fatal(err)
	}
	return dir
}

func GenerateHexString(color RGBColor) string {
	return fmt.Sprintf("%02x%02x%02x", color.R, color.G, color.B)
}
