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
	NoGap map[string]string
	AncestralNodes []string
	Length int
}

func ReadPhylip(fileName string) (alignment Alignment) {
	alignment.Alignment = make(map[string]string)
	alignment.NoGap = make(map[string]string)

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
			if strings.HasPrefix(id, "node #") {
				id = strings.Replace(id, "node #", "", -1)
				alignment.AncestralNodes = append(alignment.AncestralNodes, id)
			}
			sequence := strings.Replace(s[10:], " ", "", -1)
			alignment.NoGap[id] = strings.Replace(sequence, "-", "", -1)
			alignment.Alignment[id] = sequence
			if alignment.Length == 0 {
				alignment.Length = len(sequence)
			}
		} else {
			alignment.Header = r
		}
	}
	return alignment
}