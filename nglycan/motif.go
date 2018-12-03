package nglycan

import "regexp"

const nglycoRegex = `N[^PX][ST]`
const nRegex = `N-*`
const tsRegex = `-*[TS]`

var (
	noGapRegex = regexp.MustCompile(nglycoRegex)
	nRe        = regexp.MustCompile(nRegex)
	tsRe       = regexp.MustCompile(tsRegex)
)

// Check if a series of N and T/S position motifs do not contain P.
// Return a blank int array if no legal motifs were found.
func CheckMotif(seq string, v []int, ts [][]int) (result []int) {
	if seq[v[1]:v[1]+1] != "P" {
		for _, v1 := range ts {
			if v1[0] == v[1]+1 {
				return []int{v[0], v1[1]}
			}
		}
	}
	return result
}

// Parse N-glycosylation motif that may contain gaps within a string.
// Start by identification of all N and T/S positions. Then checking if they don't contain Proline.
// Return an [][]int containing the start and stop positions of each N-glycosylation motif.
func MotifParseStringWithGap(c string, number int) (result [][]int) {
	n := nRe.FindAllStringIndex(c, -1)
	ts := tsRe.FindAllStringIndex(c, -1)
	for _, v := range n {
		if v[1] < len(c) {
			r := CheckMotif(c, v, ts)
			if len(r) > 0 {
				result = append(result, r)
				number--
			}
			if number == 0 {
				return result
			}
		}
	}
	return result
}

func NSequoNoGap(c string) string {
	return noGapRegex.FindString(c)

}
