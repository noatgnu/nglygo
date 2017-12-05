package clustalowrapper

import (
	"github.com/noatgnu/ancestral/blastwrapper"
	"errors"
	"os/exec"
	"log"
	"os"
	"bufio"
	"io"
	"strings"
)

type ClustalOCommandline struct {
	Command string
	In string
	Out string
	OutFmt string
	blastwrapper.CommandLine
}

func (c *ClustalOCommandline) Execute() (err error) {
	commandArray, err := c.CommandBuild()
	if err != nil {
		return err
	}
	if commandArray != nil {
		cmd := exec.Command(commandArray[0], commandArray[1:]...)
		log.Printf("Started: Performing Clustal Omega Operation (%v)", c.In)
		_ = cmd.Run()
		log.Printf("Finished: Performing Clustal Omega Operation (%v)", c.In)
	}
	return err
}

func (c *ClustalOCommandline) CommandBuild() (commandArray []string, err error) {
	if c.Command != "" {
		commandArray = append(commandArray, c.Command)
	} else {
		commandArray = append(commandArray, `clustalo`)
	}
	if c.In != "" {
		commandArray = append(commandArray, "-i", c.In)
	} else {
		return nil, errors.New("text: Require Input")
	}
	if c.Out != "" {
		commandArray = append(commandArray, "-o", c.Out)
	} else {
		return nil, errors.New("text: Require Output")
	}
	if c.OutFmt != "" {
		commandArray = append(commandArray, c.OutFmt)
	}

	return commandArray, nil
}

func (c *ClustalOCommandline) ConvertSequential(fileName string, outName string) (err error) {
	f, err := os.Open(fileName)
	buff := bufio.NewReader(f)
	var header string
	count := 0
	var alignment []string
	firstPass := false
	numbSeq := 0
	enumCount := 0
	for {
		r, err := buff.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			}
			log.Panicln(err)
		}
		count++

		if count != 1 {
			s := strings.TrimSpace(r)
			if !firstPass {
				alignment = append(alignment, s)
				if s == "" {
					firstPass = true
					numbSeq = len(alignment)
				}
			}  else {
				if s == "" {
					enumCount = 0
				} else if enumCount < numbSeq {
					alignment[enumCount] += s
					enumCount++
				}
			}
		} else {
			header = r
		}
	}
	f.Close()
	f, err = os.Create(outName)
	if err != nil {
		return err
	}
	writer := bufio.NewWriter(f)
	writer.WriteString(header)
	writer.WriteString(strings.Join(alignment, "\n"))
	writer.Flush()
	f.Close()
	return err
}