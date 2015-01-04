(ns bio-algs.dynprog-test
  (:require [clojure.test :refer :all]
            [bio-algs.dynprog :refer :all]))

(deftest dyanmic-programming-functions
  (testing "does the find change function work correctly"
    (is (= 2 (find-change 24 [50 12 10 5 1])))))
