package main

import (
	"fmt"
	"math/rand"
	"os"
	"time"

	"github.com/evolbioinfo/gotree/io/newick"
	"github.com/evolbioinfo/gotree/tree"
)

func lca_main() {
	var t *tree.Tree
	var err error
	var f *os.File
	var nbtips int
	var tips []string
	var lca *tree.Node
	if f, err = os.Open(os.Args[1]); err != nil {
		panic(err)
	}
	t, err = newick.NewParser(f).Parse()
	nbtips = len(t.Tips())
	tips = []string{t.Tips()[rand.Intn(nbtips)].Name(), t.Tips()[rand.Intn(nbtips)].Name()}
	if err != nil {
		panic(err)
	}
	var start = time.Now()

	lca, _, _, _ = t.LeastCommonAncestorRooted(nil, tips...)

	fmt.Println(lca)

	elapsed := time.Since(start).Seconds()
	fmt.Printf("\nInternal time: %v\n", elapsed)

}
