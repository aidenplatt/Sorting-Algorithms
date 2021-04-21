package Sorting;

import java.util.Arrays;
import java.util.Random;

import Plotter.Plotter;

public class Sorting {

    final static int BUBBLE_VS_QUICK_LENGTH = 12;
    final static int MERGE_VS_QUICK_LENGTH = 15;
    final static int BUBBLE_VS_QUICK_SORTED_LENGTH = 12;
    final static int ARBITRARY_VS_MEDIAN_LENGTH = 16;
    final static double T = 600.0;

    /**
     * Sorts a given array using the quick sort algorithm. At each stage the pivot
     * is chosen to be the rightmost element of the subarray.
     * <p>
     * Should run in average complexity of O(nlog(n)), and worst case complexity of
     * O(n^2)
     *
     * @param arr - the array to be sorted
     */
    public static void quickSortArbitraryPivot(double[] arr) {

        quickSortArbitraryPivot(arr, 0, arr.length - 1);
    }

    /**
     * Sorts a given array using the quick sort algorithm. runs recursively, calling
     * "partition" function
     *
     * @param arr - the array to be sorted
     * @param p   - the index to begin with
     * @param r   - the index to end with
     */
    public static void quickSortArbitraryPivot(double[] arr, int p, int r) {
        if (p < r) {
            int q = partition(arr, p, r);
            quickSortArbitraryPivot(arr, p, q - 1);
            quickSortArbitraryPivot(arr, q + 1, r);
        }
    }

    /**
     * Chooses the pivot and goes through the array swapping elements as needed
     *
     * @param arr the array to be sorted
     * @param p   the index to begin with
     * @param r   the index to end with
     * @return the index of the pivot
     */
    public static int partition(double[] arr, int p, int r) {
        double pivot = arr[r];
        int i = p - 1;
        for (int j = p; j < r; j++) {
            if (arr[j] <= pivot) {
                i++;
                exchange(arr, i, j);
            }
        }
        exchange(arr, i + 1, r);
        return i + 1;
    }


    /**
     * exchanges two values of two given indexes
     *
     * @param arr the array to be sorted
     * @param i   index of first value
     * @param j   index of second value
     */
    public static void exchange(double[] arr, int i, int j) {
        double helper = arr[i];
        arr[i] = arr[j];
        arr[j] = helper;
    }

    /**
     * Sorts a given array using the quick sort algorithm. At each stage the pivot
     * is chosen in the following way: Choose 3 random elements from the array, the
     * pivot is the median of the 3 elements.
     * <p>
     * Should run in average complexity of O(nlog(n)), and worst case complexity of
     * O(n^2)
     *
     * @param arr - the array to be sorted
     */
    public static void quickSortMedianPivot(double[] arr) {
        quickSortMedianPivot(arr, 0, arr.length - 1);
    }

    /**
     * @param arr the array to be sorted
     * @param p   index of first value
     * @param r   index of last value
     */
    public static void quickSortMedianPivot(double[] arr, int p, int r) {
        if (p < r) {
            int q = Partition2(arr, p, r);
            quickSortMedianPivot(arr, p, q - 1);
            quickSortMedianPivot(arr, q + 1, r);

        }
    }

    /**
     * Chooses the pivot and goes through the array swapping elements as needed
     *
     * @param arr the array to be sorted
     * @param p   the index to begin with
     * @param r   the index to end with
     * @return the index of the pivot
     */
    public static int Partition2(double[] arr, int p, int r) {
        int pivotIndex = chooseMedian(arr, p, r);
        double pivot = arr[pivotIndex];
        arr[pivotIndex] = arr[r];
        arr[r] = pivot;
        int i = p - 1;
        for (int j = p; j < r; j++) {
            if (arr[j] <= pivot) {
                i++;
                exchange(arr, i, j);
            }
        }
        exchange(arr, i + 1, r);
        return i + 1;

    }

    /**
     * Chooses the median between three values
     *
     * @param arr the array to be sorted
     * @param p   the index to begin with
     * @param r   the index to end with
     * @return the index of the median
     */
    public static int chooseMedian(double arr[], int p, int r) {
        int mid = (p + r) / 2;
        if (arr[p] < arr[mid]) {
            if (arr[r] < arr[p]) {
                return p;
            } else if (arr[mid] < arr[r]) {
                return mid;
            } else {
                return r;
            }
        } else if (arr[r] < arr[mid]) {
            return mid;
        } else if (arr[p] < arr[r]) {
            return p;
        } else return r;

    }

    /**
     * Sorts a given array using the merge sort algorithm.
     * <p>
     * Should run in complexity O(nlog(n)) in the worst case.
     *
     * @param arr - the array to be sorted
     */

    public static void mergeSort(double[] arr) {
        mergeSort(arr, 0, arr.length - 1);

    }

    /**
     * @param arr is the arr to be sorted
     * @param p   index of first element in array
     * @param r   index of last element of array
     */
    public static void mergeSort(double[] arr, int p, int r) {
        if (p < r) {
            int q = (p + r) / 2;
            mergeSort(arr, p, q);
            mergeSort(arr, q + 1, r);
            merge(arr, p, q, r);
        }
    }

    /**
     * Divides the array to two sub arrays, and then merges them sorted
     *
     * @param arr is the array to be sorted
     * @param p   is the index to begin with
     * @param q   is the middle index of the array
     * @param r   is the index to end with
     */
    public static void merge(double[] arr, int p, int q, int r) {
        int n1 = q - p + 1;
        int n2 = r - q;
        double[] left = new double[n1 + 1];
        double[] right = new double[n2 + 1];
        for (int i = 0; i < n1; i++) {
            left[i] = arr[p + i];
        }
        for (int j = 0; j < n2; j++) {
            right[j] = arr[q + j + 1];
        }
        //Merging the sub-arrays.
        right[n2] = Double.MAX_VALUE;
        left[n1] = Double.MAX_VALUE;
        int i = 0;
        int j = 0;
        for (int k = p; k <= r; k++) {
            if (left[i] <= right[j]) {
                arr[k] = left[i];
                i++;
            } else {
                arr[k] = right[j];
                j++;
            }
        }

    }
    

    /**
     * Sorts a given array using bubble sort.
     * If at any time the algorithm recognizes no more inversions are needed it should stop.
     * <p>
     * The algorithm should run in complexity O(n^2) in the worst case.
     * The algorithm should run in complexity O(n) in the best case.
     *
     * @param arr - the array to be sorted
     */
    public static void bubbleSort(double[] arr) {
        bubbleHelper(arr, arr.length);
    }

    /**
     * Bubble helper which sorts the array using bubbleSort and returns it
     * to the main function
     *
     * @param arr is the array we sort
     * @param n   is
     */
    public static void bubbleHelper(double[] arr, int n) {
        double temp;
        int counter = 0;
        for (int i = 0; i < n; i++) {
            counter = 0;
            for (int j = 1; j < n - i; j++) {
                if (arr[j - 1] > arr[j]) {
                    temp = arr[j - 1];
                    arr[j - 1] = arr[j];
                    arr[j] = temp;
                    counter++;
                }
            }
            if ((counter == 0) && (i < n - 1)) {
//	            System.out.println("Array is ready and sorted, in linear " +
//                        "time!! :)");
                return;
            }
        }
    }

    public static void main(String[] args) {

        bubbleVsQuick();
        mergeVsQuick();
        bubbleVsQuickOnSortedArray();
        arbitraryPivotVsMedianPivot();
//		double[] arr= {2.0, 4.0, 1.0, 7.0};
//        double[] arr1= {12.0, 1.0, 7.0, 88.0, 4.0, 30.0, 29.0, 100.0};
//        double[] arr2= {4.0, 1.0};
//        double[] arr3= {4.0, 1.0, 5.0, 11.0};
//        double[] arr4 = {4.0, 4.0, 4.0, 4.0};
//        bubbleSort(arr);
//        bubbleSort(arr1);
//        bubbleSort(arr2);
//        bubbleSort(arr3);
//        bubbleSort(arr4);
    }

    /**
     * Compares the selection sort algorithm against quick sort on random arrays
     */
    public static void bubbleVsQuick() {
        double[] quickTimes = new double[BUBBLE_VS_QUICK_LENGTH];
        double[] bubbleTimes = new double[BUBBLE_VS_QUICK_LENGTH];
        long startTime, endTime;
        Random r = new Random();
        for (int i = 0; i < BUBBLE_VS_QUICK_LENGTH; i++) {
            long sumQuick = 0;
            long sumSelection = 0;
            for (int k = 0; k < T; k++) {
                int size = (int) Math.pow(2, i);
                double[] a = new double[size];
                double[] b = new double[size];
                for (int j = 0; j < a.length; j++) {
                    a[j] = r.nextGaussian() * 5000;
                    b[j] = a[j];
                }
                startTime = System.currentTimeMillis();
                quickSortArbitraryPivot(a);
                endTime = System.currentTimeMillis();
                sumQuick += endTime - startTime;
                startTime = System.currentTimeMillis();
                bubbleSort(b);
                endTime = System.currentTimeMillis();
                sumSelection += endTime - startTime;
            }
            quickTimes[i] = sumQuick / T;
            bubbleTimes[i] = sumSelection / T;
        }
        Plotter.plot("quick sort on random array", quickTimes, "bubble sort on random array", bubbleTimes);
    }

    /**
     * Compares the merge sort algorithm against quick sort on random arrays
     */
    public static void mergeVsQuick() {
        double[] quickTimes = new double[MERGE_VS_QUICK_LENGTH];
        double[] mergeTimes = new double[MERGE_VS_QUICK_LENGTH];
        long startTime, endTime;
        Random r = new Random();
        for (int i = 0; i < MERGE_VS_QUICK_LENGTH; i++) {
            long sumQuick = 0;
            long sumMerge = 0;
            for (int k = 0; k < T; k++) {
                int size = (int) Math.pow(2, i);
                double[] a = new double[size];
                double[] b = new double[size];
                for (int j = 0; j < a.length; j++) {
                    a[j] = r.nextGaussian() * 5000;
                    b[j] = a[j];
                }
                startTime = System.currentTimeMillis();
                quickSortArbitraryPivot(a);
                endTime = System.currentTimeMillis();
                sumQuick += endTime - startTime;
                startTime = System.currentTimeMillis();
                mergeSort(b);
                endTime = System.currentTimeMillis();
                sumMerge += endTime - startTime;
            }
            quickTimes[i] = sumQuick / T;
            mergeTimes[i] = sumMerge / T;
        }
        Plotter.plot("quick sort on random array", quickTimes, "merge sort on random array", mergeTimes);
    }

    /**
     * Compares the merge sort algorithm against quick sort on pre-sorted arrays
     */
    public static void bubbleVsQuickOnSortedArray() {
        double[] quickTimes = new double[BUBBLE_VS_QUICK_SORTED_LENGTH];
        double[] bubbleTimes = new double[BUBBLE_VS_QUICK_SORTED_LENGTH];
        long startTime, endTime;
        for (int i = 0; i < BUBBLE_VS_QUICK_SORTED_LENGTH; i++) {
            long sumQuick = 0;
            long sumBubble = 0;
            for (int k = 0; k < T; k++) {
                int size = (int) Math.pow(2, i);
                double[] a = new double[size];
                double[] b = new double[size];
                for (int j = 0; j < a.length; j++) {
                    a[j] = j;
                    b[j] = j;
                }
                startTime = System.currentTimeMillis();
                quickSortArbitraryPivot(a);
                endTime = System.currentTimeMillis();
                sumQuick += endTime - startTime;
                startTime = System.currentTimeMillis();
                bubbleSort(b);
                endTime = System.currentTimeMillis();
                sumBubble += endTime - startTime;
            }
            quickTimes[i] = sumQuick / T;
            bubbleTimes[i] = sumBubble / T;
        }
        Plotter.plot("quick sort on sorted array", quickTimes, "bubble sort on sorted array", bubbleTimes);
    }

    /**
     * Compares the quick sort algorithm once with a choice of an arbitrary pivot and once with a choice of a median pivot
     */
    public static void arbitraryPivotVsMedianPivot() {
        double[] arbitraryTimes = new double[ARBITRARY_VS_MEDIAN_LENGTH];
        double[] medianTimes = new double[ARBITRARY_VS_MEDIAN_LENGTH];
        long startTime, endTime;
        Random r = new Random();
        for (int i = 0; i < ARBITRARY_VS_MEDIAN_LENGTH; i++) {
            long sumArbitrary = 0;
            long sumMedian = 0;
            for (int k = 0; k < T; k++) {
                int size = (int) Math.pow(2, i);
                double[] a = new double[size];
                double[] b = new double[size];
                for (int j = 0; j < a.length; j++) {
                    a[j] = r.nextGaussian() * 5000;
                    b[j] = a[j];
                }
                startTime = System.currentTimeMillis();
                quickSortArbitraryPivot(a);
                endTime = System.currentTimeMillis();
                sumArbitrary += endTime - startTime;
                startTime = System.currentTimeMillis();
                quickSortMedianPivot(b);
                endTime = System.currentTimeMillis();
                sumMedian += endTime - startTime;
            }
            arbitraryTimes[i] = sumArbitrary / T;
            medianTimes[i] = sumMedian / T;
        }
        Plotter.plot("quick sort with an arbitrary pivot", arbitraryTimes, "quick sort with a median pivot", medianTimes);
    }

}
