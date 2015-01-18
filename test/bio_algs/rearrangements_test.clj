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
  (testing "find connected block"
    (is (= -3
           (find-connected-vertex 2 '((1 2 3 4 5 6)))
           (find-connected-vertex 2 '((1 4 5 6) (3 2)))
           (find-connected-vertex 2 '((1 -3 -2 4 5 6)))
           (find-connected-vertex -2 '((1 -3 2 4 5 6))))
        (= 3 (find-connected-vertex 2 '((1 2 -3) (4 5 6))))))
  (testing "count cycles"
    (is (= 3
           (count-cycles '((1 2) (5)) '((1 2) (5)))
           (count-cycles '((1 2 3 4 5 6)) '((1 -3 -6 -5) (2 -4))))))
  (testing "2 break distance"
    (is (= 3 (distance '((1 2 3 4 5 6)) '((1 -3 -6 -5) (2 -4))))))
  (testing "to graph and back"
    (is (= {1 -2, -2 1, 2 -3, -3 2, 3 -4, -4 3, 4 -1, -1 4} (genome-to-graph '((1 2 3 4)))))
    (is (zero? (distance '((1 2 3 4)) (graph-to-genome {1 -2, -2 1, 2 -3, -3 2, 3 -4, -4 3, 4 -1, -1 4}))))
    (is (= {1 -1, -1 1} (genome-to-graph '((1)))))
    (is (zero? (distance '((1)) (graph-to-genome {1 -1, -1 1})))))
  (testing "making a 2 break"
    (is (zero? (distance '((1 -3 -6 -5 -4 2)) (make-break '((1 -3 -6 -5) (2 -4)) '(4 -5) '(-1 2))))))
  (testing "2 break sorting"
    (let [out (two-break-sort '((1 -2 -3 4)) '((1 2 -4 -3)))]
      (is (= 4 (count out)))
      (is (= '((1 -2 -3 4)) (first out)))
      (is (= '((1 2 -4 -3)) (last out)))))
  (testing "shared kmers"
    (is (= [[0, 0] [0, 4] [4, 2] [6, 6]] (shared-kmers 3 "AAACTCATC" "TTTCAAATC")))))
