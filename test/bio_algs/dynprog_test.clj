(ns bio-algs.dynprog-test
  (:require [clojure.test :refer :all]
            [bio-algs.dynprog :refer :all]))

(deftest dyanmic-programming-functions
  (testing "does the find change function work correctly"
    (is (= 2 (find-change 24 [50 12 10 5 1])))
    (is (= 2 (find-change 40 [50 25 20 10 5 1]))))
  (testing "Manhattan Tourist Problem"
    (is (= 34
           (manhattan 4 4
                      [[1 0 2 4 3]
                       [4 6 5 2 1]
                       [4 4 5 2 1]
                       [5 6 8 5 3]]
                      [[3 2 4 0]
                       [3 2 4 2]
                       [0 7 3 3]
                       [3 3 0 2]
                       [1 3 2 2]])))))
