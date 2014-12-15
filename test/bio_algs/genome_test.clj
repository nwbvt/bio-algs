(ns bio-algs.genome-test
  (:require [clojure.test :refer :all]
            [bio-algs.genome :refer :all]) )

(deftest genome-assembly
  (testing "string spelled algorithm"
    (is (= "ACCGAAGCT" (string-spelled ["ACCGA" "CCGAA" "CGAAG" "GAAGC" "AAGCT"]))))
  (testing "overlap graph"
    (is (= #{["AGGCA", "GGCAT"], ["CATGC", "ATGCG"], ["GCATG", "CATGC"], ["GGCAT", "GCATG"]}
           (set (overlap-graph ["ATGCG" "GCATG" "CATGC" "AGGCA" "GGCAT"]))))
    (is (= ["AGGCA -> GGCAT" "CATGC -> ATGCG" "GCATG -> CATGC" "GGCAT -> GCATG"]
           (sort (format-graph (overlap-graph ["ATGCG" "GCATG" "CATGC" "AGGCA" "GGCAT"]))))))
  (testing "De Bruijn Graph"
    (is (= #{["AAG" ["AGA", "AGA"]] ["AGA" ["GAT"]] ["ATT" ["TTC"]] ["CTA" ["TAA"]] ["CTC" ["TCT"]]
             ["GAT" ["ATT"]] ["TAA" ["AAG"]] ["TCT" ["CTA" "CTC"]] ["TTC" ["TCT"]]}
           (set (deBruijn 4 "AAGATTCTCTAAGA"))))
    (is (= #{["AGG" ["GGG"]] ["CAG" ["AGG" "AGG"]] ["GAG" ["AGG"]] ["GGA" ["GAG"]] ["GGG" ["GGA" "GGG"]]}
           (set (deBruijn-kmers ["GAGG" "CAGG" "GGGG" "GGGA" "CAGG" "AGGG" "GGAG"]))))))

