package blastwrapper

import (
	"fmt"
	"regexp"
	"strings"
)

//const OrganismRegexPat = `\[([\w\s\(\)/]+)\]`
const OrganismRegexPat = `\[([^\]]+)\]`
const Species = `[\w|\s]+`
var OrgRegex = regexp.MustCompile(OrganismRegexPat)
var SpeciesRegex = regexp.MustCompile(Species)

type CommandLine interface {
	Execute() error
	CommandBuild() ([]string, error)
}

type TopDom struct {
	Start int
	Stop int
	Type string
}

type PrimeSeq struct {
	Id  string
	Seq string
	Length int
	Species string
	Name string
	TopDomain []TopDom
}

func (p *PrimeSeq) ToString() string {
	if !strings.HasPrefix(p.Id, ">") {
		return fmt.Sprintf(">%v\n%v\n", p.Id, p.Seq)
	}
	return fmt.Sprintf("%v\n%v\n", p.Id, p.Seq)
}

func SeqQualityControl(seq PrimeSeq, partialCheck bool) (quality bool) {
	if partialCheck {
		if strings.Contains(strings.ToLower(seq.Id), "partial") {
			return false
		}
	}
	if strings.Contains(seq.Id, "LOW QUALITY PROTEIN") {
		return false
	}
	for _, v := range []string{"B", "X", "Z", "J", "U"} {
		if strings.Contains(seq.Seq, v) {
			//log.Println(seq.Seq)
			return false
		}
	}
	return true
}

func SeqFilterOrganism(seq PrimeSeq, organisms []string, removeFound bool, collect map[string]bool) (check bool, organismLeft []string, found string) {
	org := ""
	if seq.Species != "" {
		org = SpeciesRegex.FindString(seq.Species)
		org = strings.TrimSpace(org)
	} else {
		org = OrgRegex.FindString(seq.Id)
	}

	//log.Println(seq.Id)
	//log.Println(org)
	if org != "" {
		if _, ok := collect[seq.Seq]; !ok {
			for _, v := range organisms {
				if strings.Contains(strings.ToLower(org), strings.ToLower(v)) {
					if removeFound {
						var arr []string
						for _, a := range organisms{
							if a != v {
								arr = append(arr, a)
							}
						}
						collect[seq.Seq] = true
						return true, arr, v
					}
					collect[seq.Seq] = true
					return true, organisms, v
				}
			}
		}

	}
	return false, nil, ""
}

func SeqFilterOrganismNoRegex(seq PrimeSeq, organisms []string) (check bool) {
	for _, v := range organisms {
		if strings.Contains(strings.ToLower(seq.Id), strings.ToLower(v)) {
			return true
		}
	}
	return false
}
