package blastwrapper

import (
	"bufio"
	"errors"
	"io"
	"log"
	"os"
	"os/exec"
	"strings"
)

type BlastDBCMDCommandline struct {
	Command        string
	In             string
	DB             string
	DBType         string
	Out            string
	FilterOrganism []string
	PartialCheck   bool
	CommandLine
}

func (b *BlastDBCMDCommandline) Execute() (err error) {
	commandArray, err := b.CommandBuild()
	if err != nil {
		return err
	}
	c := make(chan PrimeSeq)

	if commandArray != nil {
		cmd := exec.Command(commandArray[0], commandArray[1:]...)
		out, err := cmd.StdoutPipe()
		if err != nil {
			return err
		}
		go b.Process(out, c, b.FilterOrganism)
		log.Println(cmd.Args)
		err = cmd.Start()
		if err != nil {
			return err
		}
		if b.Out != "" {
			outfile, err := os.Create(b.Out)
			if err != nil {
				log.Panicln(err)
			}

			writer := bufio.NewWriter(outfile)

			for s := range c {
				writer.WriteString(s.ToString())
			}
			defer outfile.Close()
		} else {
			return errors.New("text: need output file")
		}
		defer cmd.Wait()
	} else {
		return errors.New("text: need parameters")
	}
	return err
}

func (b *BlastDBCMDCommandline) CommandBuild() (commandArray []string, err error) {
	if b.Command != "" {
		commandArray = append(commandArray, b.Command)
	} else {
		commandArray = append(commandArray, "blastdbcmd")
	}

	if b.DB != "" {
		commandArray = append(commandArray, "-db", b.DB)
	} else {
		return nil, errors.New("text: need db type")
	}
	if b.DBType != "" {
		commandArray = append(commandArray, "-dbtype", b.DBType)
	}
	commandArray = append(commandArray, "-entry_batch")
	if b.In != "" {
		commandArray = append(commandArray, b.In)
	} else {
		return nil, errors.New("text: need input file")
	}

	return commandArray, nil
}

func (b *BlastDBCMDCommandline) Process(stdout io.ReadCloser, c chan PrimeSeq, organisms []string) {
	buff := bufio.NewScanner(stdout)
	seq := PrimeSeq{"", ""}
	for buff.Scan() {
		if len(organisms) == 0 {
			break
		}
		if strings.Contains(buff.Text(), ">") {
			if seq.Id != "" {
				if SeqQualityControl(seq, b.PartialCheck) {
					check, left := SeqFilterOrganism(seq, organisms, true)
					if check {
						c <- seq
						organisms = left
					}
				}

			}
			seq = PrimeSeq{"", ""}
			seq.Id = strings.TrimSpace(buff.Text())
		} else {
			seq.Seq += strings.TrimSpace(buff.Text())
		}
	}
	if seq.Seq != "" {
		if SeqQualityControl(seq, b.PartialCheck) {
			check, left := SeqFilterOrganism(seq, organisms, true)
			if check {
				c <- seq
				organisms = left
			}
		}
	}
	close(c)
}
