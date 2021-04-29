package musclewrapper

import (
	"bufio"
	"errors"
	"fmt"
	"github.com/noatgnu/ancestral/blastwrapper"
	"io"
	"log"
	"os"
	"os/exec"
	"regexp"
	"strconv"
	"strings"
)

const phylipRow = `(\w+)\s+([\w-]+)`
var phylipRowRegex = regexp.MustCompile(phylipRow)

type MuscleCommandline struct {
	Command string
	In string
	Out string
	OutFmt string
	blastwrapper.CommandLine
}

func (c *MuscleCommandline) Execute() (err error) {
	commandArray, err := c.CommandBuild()
	if err != nil {
		return err
	}
	if commandArray != nil {
		log.Println(commandArray)
		cmd := exec.Command(commandArray[0], commandArray[1:]...)
		log.Printf("Started: Performing Clustal Omega Operation (%v)", c.In)
		_ = cmd.Run()
		log.Printf("Finished: Performing Clustal Omega Operation (%v)", c.In)
	}
	return err
}

func (c *MuscleCommandline) CommandBuild() (commandArray []string, err error) {
	if c.Command != "" {
		commandArray = append(commandArray, c.Command)
	} else {
		commandArray = append(commandArray, `clustalo`)
	}
	if c.In != "" {
		commandArray = append(commandArray, "-in", c.In)
	} else {
		return nil, errors.New("text: Require Input")
	}
	if c.Out != "" {
		if strings.HasSuffix(c.Out, "_inter"){
			commandArray = append(commandArray, "-physout", c.Out)
		} else {
			commandArray = append(commandArray, "-fastaout", c.Out)
		}

	} else {
		return nil, errors.New("text: Require Output")
	}

	return commandArray, nil
}

func (c *MuscleCommandline) ConvertAlignment(fileName string, outName string, outFormat string) (err error) {
	f, err := os.Open(fileName)
	buff := bufio.NewReader(f)
	var header string
	count := 0
	var alignment []string
	firstPass := false
	//numbSeq := 0
	enumCount := -1
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
				sSplit := strings.Split(s, " ")
				if sSplit[0] != "" {
					_, err := strconv.Atoi(sSplit[0])
					if err != nil {
						if len(s) <= 10 {
							alignment[enumCount] += strings.TrimSpace(s)
						} else {
							alignment[enumCount] += s[:10] + strings.ReplaceAll(s[10:], " ", "")
						}
					} else {
						//alignment = append(alignment, s)
						alignment = append(alignment, s[:10] + strings.ReplaceAll(s[10:], " ", ""))
						log.Println(alignment)
						enumCount++
					}
				}

				if s == "" {
					firstPass = true
					//numbSeq = len(alignment)
				}
			}
		} else {
			header = r
			//headerSplit := strings.Split(header, " ")
			//numbSeq, err = strconv.Atoi(headerSplit[0])
			//if err != nil {
			//	log.Panic(err)
			//}
			enumCount = -1
		}
	}
	f.Close()
	f, err = os.Create(outName)
	if err != nil {
		return err
	}
	writer := bufio.NewWriter(f)
	switch outFormat {
	case "phylip":
		headerSplit := strings.Split(header, " ")
		writer.WriteString(fmt.Sprintf(" %v  %v", headerSplit[0], headerSplit[1]))
		writer.WriteString(strings.Join(alignment, "\n")+"\n")

	case "fasta":
		for r := range alignment {
			a := phylipRowRegex.FindAllStringSubmatch(alignment[r], -1)
			log.Println(a)
			if len(a) >0 {
				writer.WriteString( fmt.Sprintf(">%s\n%s\n\n", a[0][1], a[0][2]))
			}
		}
	}
	writer.Flush()
	f.Close()
	return err
}