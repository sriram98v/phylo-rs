package main

import (
	"fmt"
	"os"
	"time"

	"github.com/evolbioinfo/gotree/io/newick"
	"github.com/evolbioinfo/gotree/tree"
)

func main() {
	var t *tree.Tree
	var err error
	var f *os.File
	if f, err = os.Open(os.Args[1]); err != nil {
		panic(err)
	}
	var start = time.Now()
	t, err = newick.NewParser(f).Parse()
	if err != nil {
		panic(err)
	}
	fmt.Println(t.Newick())
	elapsed := time.Since(start)
	fmt.Printf("Internal time: %v\n", elapsed)
}
