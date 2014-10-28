(ns bio-algs.core-test
  (:require [clojure.test :refer :all]
            [bio-algs.core :refer :all]))

(deftest most-common
  (testing "Most common kmer"
    (is (= '("CATG" "GCAT")
           (most-common-kmer "ACGTTGCATGTCGCATGATGCATGAGAGCT" 4)))))
