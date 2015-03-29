(ns bio-algs.patterns-test
  (:require [clojure.test :refer :all]
            [bio-algs.patterns :refer :all]))

(deftest tries
  "Testing the use of tries for pattern matching"
  (testing "Creating a trie"
     (is (= (create-trie ["ATAGA" "ATC" "GAT"])
           [[nil [1 7]]
            [\A [2]]
            [\T [3 6]]
            [\A [4]]
            [\G [5]]
            [\A []]
            [\C []]
            [\G [8]]
            [\A [9]]
            [\T []]])))
  (testing "Displaying a trie in the format the grader expects"
     (is (= (set
             (display-trie 
               [[nil [1 7]]
               ["A" [2]]
               ["T" [3 6]]
               ["A" [4]]
               ["G" [5]]
               ["A" []]
               ["C" []]
               ["G" [8]]
               ["A" [9]]
               ["T" []]]))
           #{"0->1:A"
             "1->2:T"
             "2->3:A"
             "3->4:G"
             "4->5:A"
             "2->6:C"
             "0->7:G"
             "7->8:A"
             "8->9:T"})))
  (testing "Using a trie for pattern matching"
    (let [trie (create-trie ["ATCG" "GGGT"])]
      (is (true? (match? trie "ATCGAT")))
      (is (true? (match? trie "GGGT")))
      (is (false? (match? trie "ATGC")))
      (is (= (matches trie "AATCGGGTTCAATCGGGGT")
             [1 4 11 15])))))
