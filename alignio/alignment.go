package alignio

import (
	"os"
	"log"
	"bufio"
	"io"
	"strings"
	"strconv"
)

type Alignment struct {
	Header string `json:"header"`
	Alignment map[string]string `json:"alignment"`
	NoGap map[string]string	`json:"noGap"`
	AncestralNodes []string `json:"ancestralNodes"`
	Length int `json:"length"`
	SeqIdArray []SeqD3Coordinate `json:"seqIdArray"`
	ConserveMap map[int]int `json:"conserveMap"`
}

type SeqD3Coordinate struct {
	SeqId string `json:"seqId"`
	YCoord int `json:"yCoord"`
}

type BySeqId []SeqD3Coordinate

func (s BySeqId) Len() int           { return len(s) }
func (s BySeqId) Swap(i, j int)      { s[i], s[j] = s[j], s[i] }
func (s BySeqId) Less(i, j int) bool {
	si, err := strconv.Atoi(s[i].SeqId)
	if err != nil {
		log.Panicln(err)
	}
	sj, err := strconv.Atoi(s[j].SeqId)
	if err != nil {
		log.Panicln(err)
	}
	return si < sj
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