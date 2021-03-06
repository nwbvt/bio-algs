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
    (is (= {"AGG" ["GGG"] "CAG" ["AGG" "AGG"] "GAG" ["AGG"] "GGA" ["GAG"] "GGG" ["GGG" "GGA"]}
           (deBruijn-kmers ["GAGG" "CAGG" "GGGG" "GGGA" "CAGG" "AGGG" "GGAG"]))))
  (testing "Eulerian Cycles and Paths"
    (is (= [0 3 2 6 8 7 9 6 5 4 2 1 0]
           (euler-cycle {0 [3], 1 [0], 2 [1 6], 3 [2], 4 [2], 5 [4], 6 [5 8], 7 [9], 8 [7], 9 [6]})))
    (is (= [6 3 0 2 1 3 4 6 7 8 9]
           (euler-path {0 [2], 1 [3], 2 [1], 3 [0 4], 6 [3 7], 7 [8], 8 [9], 9 [6]}))))
  (testing "String Reconstruction"
    (is (= "GGCTTACCA" (string-recon ["CTTA" "ACCA" "TACC" "GGCT" "GCTT" "TTAC"]))))
  (testing "String spelled from gapped patterns"
    (is (= "GACCGAGCGCCGGA" (gapped-string-spelled 4 2
                              [["GACC" "GCGC"]
                               ["ACCG" "CGCC"]
                               ["CCGA" "GCCG"]
                               ["CGAG" "CCGG"]
                               ["GAGC" "CGGA"]]))))
  (testing "String Reconstruction from read pairs"
    (is (= "GTGGTCGTGAGATGTTGA"
           (recon-from-read-pairs 4 2
             [["GAGA" "TTGA"]
              ["TCGT" "GATG"]
              ["CGTG" "ATGT"]
              ["TGGT" "TGAG"]
              ["GTGA" "TGTT"]
              ["GTGG" "GTGA"]
              ["TGAG" "GTTG"]
              ["GGTC" "GAGA"]
              ["GTCG" "AGAT"]])))
    (is (= "AGCAGCTGCTGCA" (recon-from-read-pairs 2 1 [["AG" "AG"]
                                                       ["GC" "GC"]
                                                       ["CA" "CT"]
                                                       ["AG" "TG"]
                                                       ["GC" "GC"]
                                                       ["CT" "CT"]
                                                       ["TG" "TG"]
                                                       ["GC" "GC"]
                                                       ["CT" "CA"]])))))
