(ns bio-algs.dynprog-test
  (:require [clojure.test :refer :all]
            [bio-algs.dynprog :refer :all]))

(deftest dyanmic-programming-functions
  (testing "does the find change function work correctly"
    (is (= 2 (find-change 24 [50 12 10 5 1])))
    (is (= 2 (find-change 40 [50 25 20 10 5 1])))))
