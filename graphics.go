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
	r, g, b byte
}

/* SoftFrameBuffer */

// SoftFrameBuffer is a simulated frame buffer that may be written to a file
type SoftFrameBuffer struct {
	Width, Height int
	buffer        [][]RGBColor
}

// NewSoftFrameBuffer allocates a new SoftFrameBuffer with specified width and
// height, with 8bit color channels
func NewSoftFrameBuffer(width, height int) *SoftFrameBuffer {
	newBuffer := &SoftFrameBuffer{
		Width:  width,
		Height: height,
		buffer: make([][]RGBColor, height, width)}

	return newBuffer
}

// WritePixel writes a pixel color to the buffer, assuming (0,0) is the
//bottom-left corner
func (frameBuffer *SoftFrameBuffer) WritePixel(x, y int, color RGBColor) {

	// We will reverse the y origin for compatibility with XPM output
	frameBuffer.buffer[frameBuffer.Height-y][x] = color

}

/* Scene */

// Scene is a virtual space that can be rendered to a frame buffer
type Scene struct {
	objects []Drawable
}

// Render renders the scene to the given buffer
func (scene *Scene) Render(buffer *SoftFrameBuffer) {
	for _, object := range scene.objects {
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
func (file *PostScriptFile) ParseObjects() ([]*Drawable, error) {

	var objects []*Drawable
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

			objects = append(objects, &object)

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
	fp, err := os.Open(path)
	if err != nil {
		return nil, err
	}

	writer := bufio.NewWriter(fp)

	file := &XPMFile{
		filePath: path,
		handle:   fp,
		writer:   writer}

	file.XPMFileFormatStr =
		"/* XPM */" +
			"static char *sco100[] = {" +
			"/* width height num_colors chars_per_pixel */" +
			"\"%d %d %d %d\"" +
			"/* colors */" +
			"%s" +
			"/* pixels */" +
			"%s" +
			"};"

	return file, err
}

func (file *XPMFile) WriteSoftBuffer(buffer *SoftFrameBuffer) {

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
	return ""
}
