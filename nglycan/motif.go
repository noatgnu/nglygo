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
