(ns bio-algs.motif-test
  (:require [clojure.test :refer :all]
            [bio-algs.motif :refer :all]))
(deftest motif-finding
  (testing "motif enumeration"
    (is (= #{"ATA" "ATT" "GTT" "TTT"}
           (motif-enum [ "ATTTGGC" "TGCCTTA" "CGGTATC" "GAAAATT" ] 3 1))))
  (testing "score function"
    (is (= 30 (score ["TCGGGGGTTTTT"
                      "CCGGTGACTTAC"
                      "ACGGGGATTTTC"
                      "TTGGGGACTTTT"
                      "AAGGGGACTTCC"
                      "TTGGGGACTTCC"
                      "TCGGGGATTCAT"
                      "TCGGGGATTCCT"
                      "TAGGGGAACTAC"
                      "TCGGGTATAACC"]))))
  (testing "distance function"
    (is (= 5 (dist "AAA" ["TTACCTTAAC"
                          "GATATCTGTC"
                          "ACGGCGTTCG"
                          "CCCTAAAGAG"
                          "CGTCAGAGGT" ]))))
  (testing "median string"
    (is (= "GAC" (median-string 3 ["AAATTGACGCAT"
                                   "GACGACCACGTT"
                                   "CGTCAGCGCCTG"
                                   "GCTGAGCACCGG"
                                   "AGTACGGGACAG"]))))
  (testing "Score kmer"
    (is (= (* 0.4 0.3 0.5 0.2 0.4)
           (score-kmer "CCGAG" {\A [0.2 0.2 0.3 0.2 0.3]
                                \C [0.4 0.3 0.1 0.5 0.1]
                                \G [0.3 0.3 0.5 0.2 0.4]
                                \T [0.1 0.2 0.1 0.1 0.2]}))))
  (testing "Making the profile"
    (is (= {\A [0.2 0.2 0.3 0.2 0.3]
            \C [0.4 0.3 0.1 0.5 0.1]
            \G [0.3 0.3 0.5 0.2 0.4]
            \T [0.1 0.2 0.1 0.1 0.2]}
           (make-profile-from-string "0.2 0.2 0.3 0.2 0.3
                                      0.4 0.3 0.1 0.5 0.1
                                      0.3 0.3 0.5 0.2 0.4
                                      0.1 0.2 0.1 0.1 0.2"))))
  (testing "generating profiles"
    (is (= {\A [0.75 1.0  0.0  0.0 ]
            \C [0.0  0.0  0.25 0.5 ]
            \G [0.0  0.0  0.0  0.25]
            \T [0.25 0.0  0.75 0.25]}
           (gen-profile ["AATC" "TATC" "AATG" "AACT"] false)))
    (is (= {\A [0.5   0.625 0.125 0.125 ]
            \C [0.125 0.125 0.25  0.375 ]
            \G [0.125 0.125 0.125 0.25  ]
            \T [0.25  0.125 0.5   0.25  ]}
           (gen-profile ["AATC" "TATC" "AATG" "AACT"] true))))
  (testing "most likely profile"
    (is (= "CCGAG" (most-probable-kmer 5
                                       "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT"
                                       {\A [0.2 0.2 0.3 0.2 0.3]
                                        \C [0.4 0.3 0.1 0.5 0.1]
                                        \G [0.3 0.3 0.5 0.2 0.4]
                                        \T [0.1 0.2 0.1 0.1 0.2]}))))
  (testing "greedy motif search"
    (is (= (score ["CAG" "CAG" "CAA" "CAA" "CAA"])
           (score (greedy-motif-search 3 ["GGCGTTCAGGCA"
                                          "AAGAATCAGTCA"
                                          "CAAGGAGTTCGC"
                                          "CACGTCAATCAC"
                                          "CAATAATATTCG"])))))
  (testing "motifs from profile function"
    (is (= ["ACCT" "ATGT" "GCGT" "ACGA" "AGGT"]
           (motifs-from-profile
             {\A [0.8 0.0 0.0 0.2]
              \C [0.0 0.6 0.2 0.0]
              \G [0.2 0.2 0.8 0.0]
              \T [0.0 0.2 0.0 0.8]}
             ["TTACCTTAAC"
              "GATGTCTGTC"
              "ACGGCGTTAG"
              "CCCTAACGAG"
              "CGTCAGAGGT"]))))
  (testing "randomized motif search"
    (is (>= (score ["TCTCGGGG" "TCTCGGGG" "TACAGGCG" "TTCAGGTG" "TCCACGTG"])
            (score (randomized-motif-search 8 ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA"
                                               "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG"
                                               "TAGTACCGAGACCGAAAGAAGTATACAGGCGT"
                                               "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC"
                                               "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"] 100))))))
(deftest entropies
  (testing "calculating entropies"
    (is (= 0.0 (entropy "AAAAAAAA")))
    (is (= 2.0 (entropy "ATCGGCTA")))
    (let [e (entropy "AAAT")]
      (is (and (< 0.811 e) (> 0.812 e)))))
  (testing "calculating motif entropies"
    (let [e (motif-entropy "ATCA" "ATGA" "ATTA" "ATAT")]
      (is (and (< 2.811 e) (> 2.812 e))))))
