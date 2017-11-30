package blastwrapper

import (
	"fmt"
	"regexp"
	"strings"
)

const OrganismRegexPat = `\[([\w\s]+)\]`

var OrgRegex = regexp.MustCompile(OrganismRegexPat)

type CommandLine interface {
	Execute() error
	CommandBuild() ([]string, error)
}

type PrimeSeq struct {
	Id  string
	Seq string
}

func (p *PrimeSeq) ToString() string {
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

	for _, v := range []string{"B", "X", "Z"} {
		if strings.Contains(seq.Seq, v) {
			return false
		}
	}
	return true
}

func SeqFilterOrganism(seq PrimeSeq, organisms []string, removeFound bool) (check bool, organismLeft []string) {
	org := OrgRegex.FindString(seq.Id)
	if org != "" {
		for i, v := range organisms {
			if strings.Contains(strings.ToLower(org), strings.ToLower(v)) {
				if removeFound {
					return true, append(organisms[:i], organisms[i+1:]...)
				} else {
					return true, organisms
				}
			}
		}
	}
	return false, nil
}

func SeqFilterOrganismNoRegex(seq PrimeSeq, organisms []string) (check bool) {
	for _, v := range organisms {
		if strings.Contains(strings.ToLower(seq.Id), strings.ToLower(v)) {
			return true
		}
	}
	return false
}
