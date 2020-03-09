/**
 * Tajonae Oyelowo
 * ID# 001888446
 */

#include <iostream>
#include <chrono>
#include <fstream>
#include <stdexcept>
const int MAX_SIZE = 1000;
using namespace std;

int moves = 0; //counts number of moves in sorting algorithms
int comps = 0; //counts number of comparisons in sorting algorithms

//Find the largest item in the array
template <class ItemType>
int findIndexofLargest(const ItemType theArray[], int size)
{
    int indexSoFar = 0; // Index of largest entry found so far
    for (int currentIndex = 1; currentIndex < size; currentIndex++)
    {
        comps++;
        // At this point, theArray[indexSoFar] >= all entries in // theArray[0..currentIndex - 1]
        if (theArray[currentIndex] > theArray[indexSoFar])
            indexSoFar = currentIndex;
    }                  //endfor
    return indexSoFar; // Index of largest entry
} // end findIndexOfLargest

//SELECTION SORT: Sorts the items in an array into ascending order.
template <class ItemType>
void selectionSort(ItemType theArray[], int n)
{
    // last = index of the last item in the subarray of items yet to be sorted;
    // largest = index of the largest item found
    for (int last = n - 1; last >= 1; last--)
    {
        // At this point, theArray[last+1..n-1] is sorted, and its entries are greater than those in theArray[0..last].
        // Select the largest entry in theArray[0..last]
        int largest = findIndexofLargest(theArray, last + 1);
        // Swap the largest entry, theArray[largest], with // theArray[last]
        if (largest != last)
        {
            swap(theArray[largest], theArray[last]);
            moves = moves + 3;
        }
    } // end for
} // end selectionSort

//BUBBLE SORT: Sorts the items in an array into ascending order.
template <class ItemType>
void bubbleSort(ItemType theArray[], int n)
{
    bool sorted = false; // False when swaps occur
    int pass = 1;
    while (!sorted && (pass < n))
    {
        // At this point, theArray[n+1-pass..n-1] is sorted
        // and all of its entries are > the entries in theArray[0..n-pass]
        sorted = true; // Assume sorted
        for (int index = 0; index < n - pass; index++)
        {
            // At this point, all entries in theArray[0..index‐1]
            // are <= theArray[index]
            int nextIndex = index + 1;
            comps++; //counts number of comparisons
            if (theArray[index] > theArray[nextIndex])
            {
                swap(theArray[index], theArray[nextIndex]);
                moves = moves + 3; //adds 3 moves to for swap
                sorted = false;    // Signal exchange
            }                      // end if
        }                          //endfor
        // Assertion: theArray[0..n‐pass‐1] < theArray[n‐pass]
        pass++;
    } // end while
} // end bubbleSort

//INSERTION SORT: Sorts the items in an array into ascending order
template <class ItemType>
void insertionSort(ItemType theArray[], int n)
{
    // unsorted = first index of the unsorted region,
    // position = index of insertion in the sorted region,
    // nextItem = next item in the unsorted region.
    // Initially, sorted region is theArray[0], unsorted region is theArray[1..n-1].
    // In general, sorted region is theArray[0..unsorted-1], unsorted region theArray[unsorted..n-1]
    for (int unsorted = 1; unsorted < n; unsorted++)
    {
        ItemType nextItem = theArray[unsorted];
        int position = unsorted;
        comps++; //counts number of comparisons
        while ((position > 0) && (theArray[position - 1] > nextItem))
        {
            // Shift theArray[position ‐ 1] to the right
            theArray[position] = theArray[position - 1];
            position--;
            comps++;           //counts number of comparisons
            moves = moves + 3; //counts number of moves
        }                      // end while
        // At this point, theArray[position] is where nextItem belongs
        if (position != unsorted)
        {                                  //to avoid unnecessary movement
            theArray[position] = nextItem; //Insert nextItem into sorted region
        }
    } //endfor
} // end insertionSort

//Merges two sorted arrays
template <class ItemType>
void merge(ItemType theArray[], int first, int mid, int last)
{
    ItemType tempArray[MAX_SIZE]; // Temporary array
    // Initialize the local indices to indicate the subarrays
    int first1 = first;   // Beginning of first subarray
    int last1 = mid;      // End of first subarray
    int first2 = mid + 1; // Beginning of second subarray
    int last2 = last;     // End of second subarray
    // While both subarrays are not empty, copy the
    //smaller item into the temporary array
    int index = first1; // Next available location in tempArray

    while ((first1 <= last1) && (first2 <= last2))
    {
        // At this point, tempArray[first..index‐1] is in order
        if (theArray[first1] <= theArray[first2])
        {
            tempArray[index] = theArray[first1];
            first1++;
            moves = moves + 3; // counts the number of moves in merge sort
        }
        else
        {
            tempArray[index] = theArray[first2];
            first2++;
            moves = moves + 3; //counts the number of moves in merge sort
        }                      // end if
        index++;
        comps++; //counts number of comparisons in merge sort
    }            // end while
    // Finish off the first subarray, if necessary
    while (first1 <= last1)
    {
        // At this point, tempArray[first..index‐1] is in order
        tempArray[index] = theArray[first1];
        first1++;
        index++;
    } // end while
    // Finish off the second subarray, if necessary
    while (first2 <= last2)
    {
        // At this point, tempArray[first..index‐1] is in order
        tempArray[index] = theArray[first2];
        first2++;
        index++;
    } //endfor
    // Copy the result back into the original array
    for (index = first; index <= last; index++)
        theArray[index] = tempArray[index];
} // end merge

//MERGE SORT: Sorts the items in an array into ascending order
template <class ItemType>
void mergeSort(ItemType theArray[], int first, int last)
{
    if (first < last)
    {
        // Sort each half
        int mid = first + (last - first) / 2; // Index of midpoint
        // Sort left half theArray[first..mid]
        mergeSort(theArray, first, mid);
        // Sort right half theArray[mid+1..last]
        mergeSort(theArray, mid + 1, last);
        // Merge the two halves
        merge(theArray, first, mid, last);
    } //endif
} // end mergeSort

//Verifies that an array is sorted is ascending order
bool sortCheck(int Array[], int size)
{
    //Checks if there is a false pair in list
    for (int i = 1; i < size; i++)
    {
        if (Array[i - 1] > Array[i])
        {
            return false;
        }
    }

    //Everything is sorted
    return true;
}

//Creates best case array
int *BSTArray()
{
    //Creates the best case array
    static int BST[1000];
    BST[0] = 10;

    //increments each element by 5 so that theyre sorted in ascending order
    for (int i = 1; i < 1000; i++)
    {
        BST[i] = BST[i - 1] + 5;
    }
    return BST;
}

//Creates average case array
int *AVGArray()
{
    //Creates the average case array
    static int AVG[1000];
    srand(time(NULL));

    //Randomizes elements of the array using rand() function
    for (int i = 0; i < 1000; i++)
    {
        AVG[i] = rand() % 100000;
    }
    return AVG;
}

//Creates worst case array
int *WSTArray()
{
    //Creates the worst case array
    static int WST[1000];
    WST[0] = 1000;

    //Decrements each element of the array so that the array is sorted in descending order
    for (int i = 1; i < 1000; i++)
    {
        WST[i] = WST[i - 1] - 1;
    }
    return WST;
}

void testSelectionSort(ofstream *ss)
{
    //Re-initializes test arrays for selection sort
    int *tBST = BSTArray();
    int *tAVG = AVGArray();
    int *tWST = WSTArray();

    //Tests best case of selection sort
    moves = 0;
    comps = 0;
    selectionSort(tBST, 1000);

    //Checks if list is sorted
    if (sortCheck(tBST, 1000))
        cout << "SS Sorted? Yes" << endl;
    else
        cout << "SS Sorted? No" << endl;

    //Writes number of moves to text file
    *ss
        << "Selection sort BST moves: " << moves << "\n";

    //Writes number of comparisons to text file
    *ss
        << "Selection sort BST comps: " << comps << "\n";

    //Tests average case of selection sort
    moves = 0;
    comps = 0;
    selectionSort(tAVG, 1000);

    //Checks if list is sorted
    if (sortCheck(tAVG, 1000))
        cout << "SS Sorted? Yes" << endl;
    else
        cout << "SS Sorted? No" << endl;

    //Writes number of moves to text file
    *ss
        << "Selection sort AVG moves: " << moves << "\n";

    //Writes number of comparisons to text file
    *ss
        << "Selection sort AVG comps: " << comps << "\n";

    //Tests worst case of selection sort
    moves = 0;
    comps = 0;
    selectionSort(tWST, 1000);

    //Checks if list is sorted
    if (sortCheck(tWST, 1000))
        cout << "SS Sorted? Yes" << endl;
    else
        cout << "SS Sorted? No" << endl;

    //Writes number of moves to text file
    *ss
        << "SS Selection sort WST moves: " << moves << "\n";

    //Writes number of comparisons to text file
    *ss
        << "SS Selection sort WST comps: " << comps << "\n";
}

void testBubbleSort(ofstream *bs)
{
    //Re-initializes test arrays for bubble sort
    int *tBST = BSTArray();
    int *tAVG = AVGArray();
    int *tWST = WSTArray();

    //Tests best case of bubble sort
    moves = 0;
    comps = 0;
    bubbleSort(tBST, 1000);

    //Checks if list is sorted
    if (sortCheck(tBST, 1000))
        cout << "BS Sorted? Yes" << endl;
    else
        cout << "BS Sorted? No" << endl;

    //Writes number of moves to text file
    *bs << "Bubble sort BST moves: " << moves << "\n";

    //Writes number of comparisons to text file
    *bs
        << "Bubble sort BST comps: " << comps << "\n";

    //Tests average case of bubble sort
    moves = 0;
    comps = 0;
    bubbleSort(tAVG, 1000);

    //Checks if list is sorted
    if (sortCheck(tAVG, 1000))
        cout << "BS Sorted? Yes" << endl;
    else
        cout << "BS Sorted? No" << endl;

    //Writes number of moves to text file
    *bs
        << "Bubble sort AVG moves: " << moves << "\n";

    //Writes number of comparisons to text file
    *bs
        << "Bubble sort AVG comps: " << comps << "\n";

    //Tests worst case of bubble sort
    moves = 0;
    comps = 0;
    bubbleSort(tWST, 1000);

    //Checks if list is sorted
    if (sortCheck(tWST, 1000))
        cout << "BS Sorted? Yes" << endl;
    else
        cout << "BS Sorted? No" << endl;

    //Writes number of moves to text file
    *bs
        << "Bubble sort WST moves: " << moves << "\n";

    //Writes number of comparisons to text file
    *bs
        << "Bubble sort WST comps: " << comps << "\n";
}

void testInsertionSort(ofstream *is)
{
    //Re-initializes test arrays for insertion sort
    int *tBST = BSTArray();
    int *tAVG = AVGArray();
    int *tWST = WSTArray();

    //Tests best case of insertion sort
    moves = 0;
    comps = 0;
    insertionSort(tBST, 1000);

    //Writes number of moves to text file
    *is
        << "Insertion sort BST moves: " << moves << "\n";

    //Writes number of comparisons to text file
    *is
        << "Insertion sort BST comps: " << comps << "\n";

    //Checks if list is sorted
    if (sortCheck(tBST, 1000))
        cout << "IS Sorted? Yes" << endl;
    else
        cout << "IS Sorted? No" << endl;

    //Tests average case of insertion sort
    moves = 0;
    comps = 0;
    insertionSort(tAVG, 1000);

    //Writes number of moves to text file
    *is
        << "Insertion sort AVG moves: " << moves << "\n";

    //Writes number of comparisons to text file
    *is
        << "Insertion sort AVG comps: " << comps << "\n";

    //Checks if list is sorted
    if (sortCheck(tAVG, 1000))
        cout << "IS Sorted? Yes" << endl;
    else
        cout << "IS Sorted? No" << endl;

    //Tests worst case of insertion sort
    moves = 0;
    comps = 0;
    insertionSort(tWST, 1000);

    //Writes number of moves to text file
    *is
        << "Insertion sort WST moves: " << moves << "\n";

    //Writes number of comparisons to text file
    *is
        << "Insertion sort WST comps: " << comps << "\n";

    //Checks if list is sorted
    if (sortCheck(tWST, 1000))
        cout << "IS Sorted? Yes" << endl;
    else
        cout << "IS Sorted? No" << endl;
}

void testMergeSort(ofstream *ms)
{
    //Re-initializes test arrays for merge sort
    int *tBST = BSTArray();
    int *tAVG = AVGArray();
    int *tWST = WSTArray();

    //Tests best case of merge sort
    moves = 0;
    comps = 0;
    mergeSort(tBST, 0, MAX_SIZE - 1);

    //Writes number of moves to text file
    *ms
        << "Merge sort BST moves: " << moves << "\n";

    //Writes number of comparisons to text file
    *ms
        << "Merge sort BST comps: " << comps << "\n";

    //Checks if list is sorted
    if (sortCheck(tBST, 1000))
        cout << "MS Sorted? Yes" << endl;
    else
        cout << "MS Sorted? No" << endl;

    //Tests average case of merge sort
    moves = 0;
    comps = 0;
    mergeSort(tAVG, 0, MAX_SIZE - 1);

    //Writes number of moves to text file
    *ms
        << "Merge sort AVG moves: " << moves << "\n";

    //Writes number of comparisons to text file
    *ms
        << "Merge sort AVG comps: " << comps << "\n";

    //Checks if list is sorted
    if (sortCheck(tAVG, 1000))
        cout << "MS Sorted? Yes" << endl;
    else
        cout << "MS Sorted? No" << endl;

    //Tests worst case of merge sort
    moves = 0;
    comps = 0;
    mergeSort(tWST, 0, MAX_SIZE - 1);

    //Writes number of moves to text file
    *ms
        << "Merge sort WST moves: " << moves << "\n";

    //Writes number of comparisons to text file
    *ms
        << "Merge sort WST comps: " << comps << "\n";

    //Checks if list is sorted
    if (sortCheck(tWST, 1000))
        cout << "MS Sorted? Yes" << endl;
    else
        cout << "MS Sorted? No" << endl;
}

int main()
{
    //Creates sort.txt file
    ofstream outf;
    outf.open("sort.txt");

    if (outf.fail())
    {
        cerr << "Error: Could not open output file\n";
        exit(1);
    }

    //Runs selection sort for three arrays and writes moves to sort.txt file
    testSelectionSort(&outf);
    testBubbleSort(&outf);
    testInsertionSort(&outf);
    testMergeSort(&outf);
    outf.close(); //Close the file at the end of your program.

    return 0;
}
