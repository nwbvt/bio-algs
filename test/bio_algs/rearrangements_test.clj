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
    (is (= 8 (count-breakpoints [3 4 5 -12 -8 -7 -6 1 2 10 9 -11 13 14]))))
  (testing "find next block"
    (is (= 3
           (find-next-block 2 '((1 2 3 4 5 6)))
           (find-next-block 2 '((1 4 5 6) (3 2)))
           (find-next-block 2 '((1 3 -2 4 5 6))))
        (= -3 (find-next-block 2 '((1 2 -3) (4 5 6))))))
  #_(testing "2 break distance"
    (is (= 3 (two-break '((1 2 3 4 5 6)) '((1 -3 -6 -5) (2 -4)))))))
