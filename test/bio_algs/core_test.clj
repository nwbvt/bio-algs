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
