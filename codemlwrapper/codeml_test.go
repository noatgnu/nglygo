package codemlwrapper

import (
	"testing"
	"log"
	"bufio"
	"os"
	"io"
	"strings"
)

func TestCodeMLCommandline_ReadSupplemental(t *testing.T) {

	f, err := os.Open("./rst")
	if err != nil {
		log.Panicln(err)
	}
	buff := bufio.NewReader(f)
	bf, err := os.Create("./rst"+"_branch")
	if err != nil {
		log.Panicln(err)
	}

	writer := bufio.NewWriter(bf)
	writer.WriteString("Origin\tTarget\tPosition\tOrigin_Residue\tStats\tTarget_Residue\n")
	//branches = make(map[string]*Branch)
	lineChan := make(chan string)
	go func() {
		for {
			r, err := buff.ReadString('\n')
			if err == io.EOF {
				break
			}
			if err != nil {
				log.Panicln(err)
			}
			s := strings.TrimSpace(string(r))
			s = strings.TrimSpace(s)
			lineChan <- s
		}
		close(lineChan)
	}()
	for s := range lineChan {
		if s != "" {
			if strings.HasPrefix(s,"tree with node labels") {
				treeFileName := "./rst_tree"
				f, err := os.Create(treeFileName)
				if err != nil {
					log.Panicln(err)
				}
				writer := bufio.NewWriter(f)
				tree := ReadNextLine(lineChan)
				writer.WriteString(tree+"\n")
				writer.Flush()
				f.Close()
				log.Printf("Wrote Reconstructed Tree %v", treeFileName)

			} else if strings.HasPrefix(s, "List of extant and reconstructed sequences") {
				log.Printf("Getting Reconstructed Alignment %v", "./rst_tree")
			}
		}
	}

	writer.Flush()
	f.Close()
	bf.Close()
}
