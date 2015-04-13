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

(deftest suffix-trees
  (testing "creating a suffix tree"
    (is (= (suffix-tree "ATAAATG$")
           [[nil [1 2 11 12]]
            [[0 1] [3 4]]
            [[1 1] [9 10]]
            [[1 1] [7 8]]
            [[3 1] [5 6]]
            [[4 4] []]
            [[5 3] []]
            [[2 6] []]
            [[6 2] []]
            [[2 6] []]
            [[6 2] []]
            [[6 2] []]
            [[7 1] []]])))
  (testing "displaying the labels of a suffix tree"
    (is (= (sort 
             (labels "ATAAATG$"
               [[nil [1 2 11 12]]
                [[0 1] [3 4]]
                [[1 1] [9 10]]
                [[1 1] [7 8]]
                [[3 1] [5 6]]
                [[4 4] []]
                [[5 3] []]
                [[2 6] []]
                [[6 2] []]
                [[2 6] []]
                [[6 2] []]
                [[6 2] []]
                [[7 1] []]]))
           (sort
             ["AAATG$"
              "G$"
              "T"
              "ATG$"
              "TG$"
              "A"
              "A"
              "AAATG$"
              "G$"
              "T"
              "G$"
              "$"]))))
  (testing "using suffix trees for various problems"
    (is (= (longest-repeat "ATATCGTTTTATCGTT")
           "TATCGTT"))
    (is (= (longest-common "TCGGTAGATTGCGCCCACTC" "AGGGGCTCGCAGTGTAAGAA") "AGA"))
    (is (true? (match-2? "ATAAATG$" "AAA")))
    (is (true? (match-2? "ATAAATG$" "TG$")))
    (is (false? (match-2? "ATAAATG$" "AAAA")))
    (is (false? (match-2? "ATAAATG$" "AATA")))
    (is (= (shortest-unique "CCAAGCTGCTAGAGG$" "CATGCTGGGCTGGCT$") "CC"))))

(deftest suffix-arrays
  (testing "Creating a suffix array"
    (is (= (suffix-array "AACGATAGCGGTAGA$")
           [15, 14, 0, 1, 12, 6, 4, 2, 8, 13, 3, 7, 9, 10, 11, 5]))))

(deftest burrows-wheeler-transforms
  (testing "Creating a Burrows-Wheeler Transform"
    (is (= (bw-transform "GCGTGCCTGGTCA$")
           "ACTGGCT$TGCGGC")))
  (testing "Reforming from a Burrows-Wheeler Transform"
    (is (= (bw-recon "TTCCTAACG$A")
           "TACATCACGT$")))
  (testing "BW matching counts"
    (let [bw "TCCTCTATGAGATCCTATTCTATGAAACCTTCA$GACCAAAATTCTCCGGC"]
      (is (= (map (partial bw-match-count bw)["CCT" "CAC" "GAG" "CAG" "ATC"])
             [2 1 1 0 1]))))
  (testing "BW matching with positions"
    (let [text "AATCGGGTTCAATCGGGGT$"]
      (is (= (sort (bw-match text ["ATCG" "GGGT"]))
             [1 4 11 15])))))
