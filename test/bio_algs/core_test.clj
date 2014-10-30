(ns bio-algs.core-test
  (:require [clojure.test :refer :all]
            [bio-algs.core :refer :all]))

(deftest most-common
  (testing "Most common kmer"
    (is (= (set ["CATG" "GCAT"])
           (set (most-common-kmer "ACGTTGCATGTCGCATGATGCATGAGAGCT" 4))))))

(deftest complement
  (testing "Finding the reverse complement"
    (is (= "ACCGGGTTTT"
           (reverse-comp "AAAACCCGGT")))))

(deftest pattern-match
  (testing "Finding a pattern"
    (is (= [1 3 9]
           (pat-match "ATAT" "GATATATGCATATACTT")))))

(deftest finding-clumps
  (testing "Finding clumps"
    (is (= #{"CGACA" "GAAGA"}
           (find-clumps "CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA"
                        5 50 4))))
  (testing "Boundary condition"
    (is (= #{} (find-clumps "ATCCATC" 3 6 2)))
    (is (= #{"ATC"} (find-clumps "ATCCATCGA" 3 7 2)))))

(deftest skew-counts
  (testing "Skew counts"
    (is (= [0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2]
          (skew "CATGGGCATCGGCCATACGCC")))))

(deftest minimum-skew
  (testing "Finding the minimum skew positions"
    (is (= [11 24]
           (min-skew "TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT")))))

(deftest hamming-distance
  (testing "Finding the hamming distance"
    (is (= 3 (hamming-dist "GGGCCGTTGGT" "GGACCGTTGAC")))))
