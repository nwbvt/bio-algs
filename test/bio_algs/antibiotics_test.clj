(ns bio-algs.antibiotics-test
  (:require [clojure.test :refer :all]
            [bio-algs.antibiotics :refer :all]))

(deftest protein-translation
  (testing "That dna translates into proteins correctly"
    (is (= "MAMAPRTEINSTRING"
           (translate "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA")))))

(deftest peptide-encoding
  (testing "Finding encodings for a peptide"
    (is (= ["ATGGCC" "GGCCAT" "ATGGCC"]
           (find-encoding "ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA" "MA")))))

(deftest theoretical-spectrum-test
  (testing "Finding the theoretical spectrum for a peptide"
    (is (= [0 113 114 128 129 227 242 242 257 355 356 370 371 484]
           (theoretical-spectrum "LEQN")))))

(deftest branch-and-bounds-test
  (testing "Using branch and bounds to perform sequencing"
    (is (= ["186-128-113" "186-113-128" "128-186-113" "128-113-186" "113-186-128" "113-128-186" ]
          (format-bb-seq [0 113 128 186 241 299 314 427]) ))))


(deftest cyclopeptide-scoring-test
  (testing "Scoring a peptide's experimental spectrum"
    (is (= 11 (score "NQEL" [0 99 113 114 128 227 257 299 355 356 370 371 484])))))
