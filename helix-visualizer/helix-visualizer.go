package main

import (
    "bufio"
    "fmt"
    "io"
    "os"
    "strings"
    "time"
)

// Nucleotide colors
const (
    Reset   = "\033[0m"
    Black   = "\033[40m"
    Blue    = "\033[44m"
    Purple  = "\033[45m"
    Red     = "\033[41m"
    Green   = "\033[42m"
    Blink   = "\033[5m"
    Grey    = "\033[90m"
    GreyEnd = "\033[0m"
)

// Nucleotides mapping
var Nucleotides = map[byte]string{
    'A': Green + "A" + Reset,
    'T': Red + "T" + Reset,
    'C': Blue + "C" + Reset,
    'G': Black + "G" + Reset,
    'N': Red + "N" + Reset,
}

func getStructure(n byte, i int, colorNuc bool, colorRevComp bool) string {
    r := complement(n)

    nucStructure := []string{
        "       %s%s   ",
        "     %s--%s   ",
        "  %s====---%s ",
        " %s=====---%s ",
        " %s---=====%s ",
        "  %s--====%s  ",
        "    %s--%s    ",
        "     %s%s     ",
        "    %s%s      ",
        "    %s--%s    ",
        "  %s--====%s  ",
        " %s---=====%s ",
        " %s=====---%s ",
        "  %s====--%s  ",
        "    %s--%s    ",
        "      %s%s    ",
    }

    nucStr := fmt.Sprintf(nucStructure[i%16], colorNucOrPlain(n, colorNuc), colorNucOrPlain(r, colorRevComp))
    nucStr = strings.ReplaceAll(nucStr, "-", Grey+"-"+Reset)
    nucStr = strings.ReplaceAll(nucStr, "=", Grey+"="+Reset)

    return nucStr
}

func colorNucOrPlain(n byte, colorNuc bool) string {
    if colorNuc {
        return Nucleotides[n]
    }
    return string(n)
}

func complement(n byte) byte {
    switch n {
    case 'A':
        return 'T'
    case 'T':
        return 'A'
    case 'C':
        return 'G'
    case 'G':
        return 'C'
    default:
        return n
    }
}

func main() {
    if len(os.Args) != 2 {
        fmt.Printf("Usage: %s <GENBANK file>\n", os.Args[0])
        os.Exit(1)
    }

    inputFile := os.Args[1]

    file, err := os.Open(inputFile)
    if err != nil {
        fmt.Printf("Error opening file: %v\n", err)
        os.Exit(1)
    }
    defer file.Close()

    reader := bufio.NewReader(file)

    for {
        line, _, err := reader.ReadLine()
        if err == io.EOF {
            break
        }
        if err != nil {
            fmt.Printf("Error reading file: %v\n", err)
            os.Exit(1)
        }

    if strings.HasPrefix(string(line), ">") {
            // This is a FASTA header, you can customize how you handle it.
            fmt.Println(string(line))
        } else {
            // This is a sequence line, process it as before.
            for i, nucleotide := range line {
                fmt.Printf("%s\n", getStructure(byte(nucleotide), i, true, true))
                time.Sleep(50 * time.Millisecond)
            }
            //fmt.Println()

            //time.Sleep(100 * time.Millisecond)
            // fmt.Print("\033[H\033[2J") // Clear the screen
        }

    }
}
