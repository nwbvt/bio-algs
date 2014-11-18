(ns bio-algs.motif-test
  (:require [clojure.test :refer :all]
            [bio-algs.motif :refer :all]))

(deftest motif-finding
  (testing "motif enumeration"
    (is (= ["ATA" "ATT" "GTT" "TTT"]
           (motif-enum {"ATTTGGC" "TGCCTTA" "CGGTATC" "GAAAATT"} 3 1))))) 
