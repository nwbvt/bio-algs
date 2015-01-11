(ns bio-algs.dynprog-test
  (:require [clojure.test :refer :all] [bio-algs.dynprog :refer :all])) 
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
                       [1 3 2 2]]))))
  (testing "Longest Common Subsequence"
    (is (= [6 "AAC-CT-TGG" "-ACACTGTGA"] (longest-common-subseq "AACCTTGG" "ACACTGTGA")))
    (is (= [8 "PLEASANTLY" "-ME--AN-LY"] (longest-common-subseq "PLEASANTLY" "MEANLY" 5 blosum62))) 
    (is (= [15 "EANL-Y" "ENALTY"] (longest-common-subseq "MEANLYM" "PENALTYP" 5 pam250 :local)))  
    (is (= [2 "TAGGCTTA" "TA-G-ATA"] (longest-common-subseq "GTAGGCTTAAGGTTA" "TAGATA" 1 (partial simple-matching 1) :fit))))   
  (testing "Topological ordering" (is (= [0 3 4 1 2]
           (top-order [[0 1]
                       [1 2]
                       [3 1]
                       [3 4]
                       [4 1]
                       [4 2]]))))
  (testing "Longest Common Path"
    (is (= [9 [0 2 3 4]]
           (longest-common-path
             0 4
             [[0 1 7]
              [0 2 4]
              [2 3 2]
              [1 4 1]
              [3 4 3]]))))
  (testing "Edit distance"
    (is (= 5
           (edit-distance "PLEASANTLY" "MEANLY")
           (edit-distance "MEANLY" "PLEASANTLY")))))
