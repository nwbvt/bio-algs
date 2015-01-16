(ns bio-algs.rearrangements-test
  (:require [clojure.test :refer :all]
            [bio-algs.rearrangements :refer :all]))

(deftest rearrangements
  (testing "reversals"
    (= [1 2 3 4 5 6] (reversal [1 2 3 -5 -4 6] 3 5)))
  (testing "greedy sorting"
    (is (= [[-1 -4 3 5 -2] [1 -4 3 5 -2] [1 2 -5 -3 4] [1 2 3 5 4] [1 2 3 -4 -5] [1 2 3 4 -5] [1 2 3 4 5]]
           (greedy-sorting [-3 4 1 5 -2]))))
  (testing "formatting permutations"
    (is (= "(-1 -4 +3 +5 -2)" (format-perm [-1 -4 3 5 -2]) (format-perm '(-1 -4 3 5 -2)))))
  (testing "counting breakpoints"
    (is (= 8 (count-breakpoints [3 4 5 -12 -8 -7 -6 1 2 10 9 -11 13 14])))))
