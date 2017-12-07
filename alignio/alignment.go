package alignio

import (
	"os"
	"log"
	"bufio"
	"io"
	"strings"
)

type Alignment struct {
	Header string
	Alignment map[string]string
}

func ReadPhylip(fileName string) (alignment Alignment) {
	alignment.Alignment = make(map[string]string)
	f, err := os.Open(fileName)
	if err != nil {
		log.Panicln(err)
	}
	defer f.Close()
	reader := bufio.NewReader(f)
	for {
		r, err := reader.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			}
			log.Panicln(err)
		}
		if alignment.Header != "" {
			s := strings.TrimSpace(r)
			id := strings.TrimSpace(s[0:10])
			sequence := strings.Replace(s[10:], " ", "", -1)
			alignment.Alignment[id] = sequence
		} else {
			alignment.Header = r
		}
	}
	return alignment
}