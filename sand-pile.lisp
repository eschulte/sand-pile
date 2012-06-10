;;; sand-pile.lisp --- 2D CA simulating a sand pile

;; Copyright 2012 Eric Schulte

;; This program is free software: you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation, either version 3 of the License, or
;; (at your option) any later version.

;; This program is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.

;; You should have received a copy of the GNU General Public License
;; along with this program.  If not, see <http://www.gnu.org/licenses/>.

;;; Commentary:

;; This program implements the 2D CA used to investigate
;; self-organized criticality through the simulation of a simplified
;; sand-pile model in the following publication.

;; @article{bak1988self,
;;   author={Bak, P. and Tang, C. and Wiesenfeld, K. and others},
;;   title={Self-organized criticality},
;;   journal={Physical review A},
;;   year={1988},
;;   volume={38},
;;   number={1},
;;   pages={364--374}
;; }

;;; Code:
(require 'alexandria)
(defpackage #:sand-pile (:use :common-lisp :alexandria))
(in-package :sand-pile)

(defvar max-init       100)
(defvar critical-slope 2)
(defvar *size*         40)
(defvar *dim*          2)
(defvar *grid*
  (make-array (make-list *dim* :initial-element *size*) :element-type :number))

(defun z (coord) (apply #'aref *grid* coord))
(defun set-z (coord new) (setf (apply #'aref *grid* coord) new))
(defsetf z set-z)

(defun range (n)
  (loop for i upto (1- n) collect i))

(defun cross (dimensions)
  (if (cdr dimensions)
      (let ((this (range (car dimensions))))
        (mapcan (lambda (rst) (mapcar (lambda (el) (cons el rst)) this))
                (cross (cdr dimensions))))
      (mapcar #'list (range (car dimensions)))))

(defun neighbors (dim)
  (remove-if (lambda (dim) (some (lambda (n) (or (< n 0) (<= *size* n))) dim))
             (mapcar (lambda (off) (map 'list #'+ off dim))
                     (mapcar (lambda (off) (mapcar #'1- off))
                             (remove (make-list *dim* :initial-element 1)
                               (cross (make-list *dim* :initial-element 3))
                               :test #'tree-equal)))))

(defun init ()
  "Randomly initialize the `*grid*'."
  (mapc (lambda (coord) (setf (z coord) (random max-init)))
        (cross (array-dimensions *grid*))))

(defun update (coord)
  "Update a coordinate in the `*grid*'."
  (flet ((slide (a b) (when (> (- (z a) (z b)) critical-slope)
                        (incf (z b)) (decf (z a)) t)))
    (mapc #'update
          (remove-if-not (lambda (neighbor)
                           (when (or (slide coord neighbor)
                                     (slide neighbor coord))
                             neighbor))
                         (sort (neighbors coord) #'<
                               :key (lambda (it) (apply #'aref *grid* it)))))))

(defun gnuplot (stream title)
  (unless (= *dim* 2) (error "can only plot 2D grids"))
  ;; Not used in interactive runs, used to generate the movie below.
  ;; (format stream "set term png~%")
  ;; (format stream "set output '~~/src/sand-pile/plots/~a.png'~%" title)
  (format stream "set title '~a'~%" title)
  (format stream "splot '-' matrix with pm3d notitle~%")
  (format stream "~{~{~a~^ ~}~%~}e~%EOF~%"
          (loop :for x :below (array-dimension *grid* 0)
             :collect (loop :for y :below (array-dimension *grid* 1)
                         :collect (aref *grid* x y)))))

(defun run (steps stream)
  (init)
  (dotimes (n steps)
    (let ((coord (random-elt (cross (array-dimensions *grid*)))))
      (incf (z coord)) (update coord))
    (gnuplot stream n)))
