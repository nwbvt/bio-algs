(ns bio-algs.phylogeny-test
  (:require [clojure.test :refer :all]
            [bio-algs.phylogeny :refer :all]))

(deftest weights
  (testing "parsing the graph"
    (= (parse-graph ["0->1:2" "0->2:3" "1->0:2" "1->3:4" "2->0:3" "3->1:4"])
       {:0 {:1 2 :2 3}
        :1 {:0 2 :3 4}
        :2 {:0 3}
        :3 {:1 4}}))
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
           [[0	13	21	22]
            [13	0	12	13]
            [21	12	0	13]
            [22	13	13	0]]))))
