(ns bio-algs.phylogeny-test
  (:require [clojure.test :refer :all]
            [bio-algs.phylogeny :refer :all]))

(deftest weights
  (testing "parsing the graph"
    (is (= (parse-graph ["0->1:2" "0->2:3" "1->0:2" "1->3:4" "2->0:3" "3->1:4"])
           {0 {1 2.0 2 3.0}
            1 {0 2.0 3 4.0}
            2 {0 3.0}
            3 {1 4.0}})))
  (testing "find the weight matrix"
    (is (= (weight-matrix (parse-graph ["0->4:11"
                                        "1->4:2"
                                        "2->5:6"
                                        "3->5:7"
                                        "4->0:11"
                                        "4->1:2"
                                        "4->5:4"
                                        "5->4:4"
                                        "5->3:7"
                                        "5->2:6"]))
           [[0.0 13.0 21.0 22.0]
            [13.0 0.0 12.0 13.0]
            [21.0 12.0 0.0 13.0]
            [22.0 13.0 13.0 0.0]]))))

(deftest finding-limb-length
  (testing "Finding the limb length of a node to its parent from a weight matrix"
    (is (= (limb-length [[0  13 21 22]
                         [13 0  12 13]
                         [21 12 0  13]
                         [22 13 13 0 ]] 1)
           2))))

(deftest adaptive-phylogeny-algorithm
  (testing "Creating the tree matching the given distance matrix"
    (is (= (make-tree [[0.0 13.0 21.0 22.0]
                       [13.0 0.0 12.0 13.0]
                       [21.0 12.0 0.0 13.0]
                       [22.0 13.0 13.0 0.0]])
           (parse-graph ["0->4:11"
                         "1->4:2"
                         "2->5:6"
                         "3->5:7"
                         "4->0:11"
                         "4->1:2"
                         "4->5:4"
                         "5->4:4"
                         "5->3:7"
                         "5->2:6"])))))

(deftest upgma-test
  (testing "Implementation of UPGMA algoritm"
    (is (= (upgma [[0	20	17	11]
                   [20	0	20	13]
                   [17	20	0	10]
                   [11	13	10	0 ]])
           (parse-graph ["0->5:7.000"
                         "1->6:8.833333333333332"
                         "2->4:5.000"
                         "3->4:5.000"
                         "4->2:5.000"
                         "4->3:5.000"
                         "4->5:2.000"
                         "5->0:7.000"
                         "5->4:2.000"
                         "5->6:1.833333333333333"
                         "6->5:1.833333333333333"
                         "6->1:8.833333333333332"])))))

(deftest neighbor-joining
  (testing "D* generation"
    (is (= (d* [[0  13 21 22]
                [13 0  12 13]
                [21 12 0  13]
                [22 13 13 0 ]])
           [[0   -68 -60 -60]
            [-68 0   -60 -60]
            [-60 -60 0   -68]
            [-60 -60 -68 0  ]])))
  (testing "Implmentation of the neighbor joining algoirthm"
    (is (= (nj-tree [[0	 23 27 20]
                     [23 0  30 28]
                     [27 30 0  30]
                     [20 28 30 0]])
           (parse-graph ["0->4:8.000"
                         "1->5:13.500"
                         "2->5:16.500"
                         "3->4:12.000"
                         "4->5:2.000"
                         "4->0:8.000"
                         "4->3:12.000"
                         "5->1:13.500"
                         "5->2:16.500"
                         "5->4:2.000"])))))

(deftest small-parsimony-problem
  (testing "Implementation of the small parsimony algorithm"
    (let [tree (parse-graph ["4->CAAATCCC"
                             "4->ATTGCGAC"
                             "5->CTGCGCTG"
                             "5->ATGGACGA"
                             "6->4"
                             "6->5"])
          mapping (small-parsimony-mapping tree)]
     (is (= mapping {4 "ATAGCCAC" 5 "ATGGACGA" 6 "ATAGACAA"}))
     #_(is (= (total-cost graph mapping) 16))
     #_(is (= 
           [16 (parse-graph ["ATTGCGAC->ATAGCCAC:2"
                             "ATAGACAA->ATAGCCAC:2"
                             "ATAGACAA->ATGGACTA:2"
                             "ATGGACGA->ATGGACTA:1"
                             "CTGCGCTG->ATGGACTA:4"
                             "ATGGACTA->CTGCGCTG:4"
                             "ATGGACTA->ATGGACGA:1"
                             "ATGGACTA->ATAGACAA:2"
                             "ATAGCCAC->CAAATCCC:5"
                             "ATAGCCAC->ATTGCGAC:2"
                             "ATAGCCAC->ATAGACAA:2"
                             "CAAATCCC->ATAGCCAC:5"])])))))
