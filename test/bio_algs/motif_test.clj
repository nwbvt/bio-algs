(ns bio-algs.motif-test
  (:require [clojure.test :refer :all]
            [bio-algs.motif :refer :all]))

(deftest motif-finding
  (testing "motif enumeration"
    (is (= #{"ATA" "ATT" "GTT" "TTT"}
           (motif-enum [ "ATTTGGC" "TGCCTTA" "CGGTATC" "GAAAATT" ] 3 1))))) 

(deftest entropies
  (testing "calculating entropies"
    (is (= 0.0 (entropy "AAAAAAAA")))
    (is (= 2.0 (entropy "ATCGGCTA")))
    (let [e (entropy "AAAT")]
      (is (and (< 0.811 e) (> 0.812 e)))))
  (testing "calculating motif entropies"
    (let [e (motif-entropy "ATCA" "ATGA" "ATTA" "ATAT")]
      (is (and (< 2.811 e) (> 2.812 e))))))

(deftest medians
  (testing "distance function"
    (is (= 5 (dist "AAA" ["TTACCTTAAC"
                          "GATATCTGTC"
                          "ACGGCGTTCG"
                          "CCCTAAAGAG"
                          "CGTCAGAGGT" ])))
    )
  (testing "median string"
    (is (= "GAC" (median-string 3 ["AAATTGACGCAT"
                                   "GACGACCACGTT"
                                   "CGTCAGCGCCTG"
                                   "GCTGAGCACCGG"
                                   "AGTACGGGACAG"])))))
