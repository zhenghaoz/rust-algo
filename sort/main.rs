extern crate rand;

use rand::Rng;
use std::time::Instant;
use std::collections::VecDeque;

struct SortAlgorithm {
    name: String,
    func: fn(&mut Vec<usize>)
}

fn bubble_sort(nums: &mut Vec<usize>) {
    for _i in 1..nums.len() {
        for i in 1..nums.len() {
            if nums[i-1] > nums[i] {
                nums.swap(i-1, i);
            }
        }
    }
}

fn insert_sort(nums: &mut Vec<usize>) {
    for i in 1..nums.len() {
        let (mut p, v) = (i, nums[i]);
        while p > 0 && nums[p-1] > v {
            nums[p] = nums[p-1];
            p -= 1;
        }
        nums[p] = v;
    }
}

fn selection_sort(nums: &mut Vec<usize>) {
    for i in 0..nums.len()-1 {
        let mut p = i;
        for j in i+1..nums.len() {
            if nums[j] < nums[p] {
                p = j;
            }
        }
        nums.swap(i, p);
    }
}

fn shell_sort(nums: &mut Vec<usize>) {

    fn _insert_sort(nums: &mut Vec<usize>, g: usize) {
        for i in g..nums.len() {
            let (mut p, v) = (i, nums[i]);
            while p >= g && nums[p-g] > v {
                nums[p] = nums[p-g];
                p -= g;
            }
            nums[p] = v;
        }
    }

    let mut a: VecDeque<usize> = VecDeque::new();
    a.push_front(1);
    while *a.front().unwrap() <= nums.len() {
        a.push_front(3*a.front().unwrap()+1);
    }
    for &g in a.iter() {
        _insert_sort(nums, g);
    }
}

fn quick_sort(nums: &mut Vec<usize>) {

    fn _partition(nums: &mut Vec<usize>, begin: usize, end: usize) -> usize {
        let (mut i, v) = (begin, nums[end-1]);
        for j in begin..end-1 {
            if nums[j] <= v {
                nums.swap(i, j);
                i += 1;
            }
        }
        nums.swap(i, end-1);
        i
    }

    fn _quick_sort(nums: &mut Vec<usize>, begin: usize, end: usize) {
        if begin+1 < end {
            let mid = _partition(nums, begin, end);
            _quick_sort(nums, begin, mid);
            _quick_sort(nums, mid+1, end);
        }
    }

    _quick_sort(nums, 0, nums.len())
}

fn merge_sort(nums: &mut Vec<usize>) {

    fn _merge(nums: &mut Vec<usize>, left: usize, mid: usize, right: usize) {
        let left_part: Vec<usize> = nums[left..mid].iter().cloned().collect();
        let right_part: Vec<usize> = nums[mid..right].iter().cloned().collect();
        let (mut left_offset, mut right_offset) = (0usize, 0usize);
        while left_offset < left_part.len() || right_offset < right_part.len() {
            if right_offset == right_part.len() 
            || (left_offset < left_part.len() && left_part[left_offset] <= right_part[right_offset]) {
                nums[left + left_offset + right_offset] = left_part[left_offset];
                left_offset += 1;
            } else {
                nums[left + left_offset + right_offset] = right_part[right_offset];
                right_offset += 1;
            }
        }
    }

    fn _merge_sort(nums: &mut Vec<usize>, left: usize, right: usize) {
        if left+1 < right {
            let mid = (left + right) / 2;
            _merge_sort(nums, left, mid);
            _merge_sort(nums, mid, right);
            _merge(nums, left, mid, right);
        }
    }

    _merge_sort(nums, 0, nums.len())
}

fn count_sort(nums: &mut Vec<usize>) {
    let n = nums.iter().max().unwrap();
    let origin_nums = nums.clone();
    let mut count: Vec<usize> = Vec::new();
    for _i in 0..n+1 {
        count.push(0)
    }
    for &v in nums.iter() {
        count[v] += 1;
    }
    for i in 1..count.len() {
        count[i] += count[i-1];
    }
    for &v in origin_nums.iter() {
        nums[count[v]-1] = v;
        count[v] -= 1;
    }
}

fn is_sorted(nums: &Vec<usize>) -> bool {
    for i in 1..nums.len() {
        if nums[i-1] > nums[i] {
            return false;
        }
    }
    true
}

struct Heap<T: Ord> {
    elems: Vec<T>
}

impl<T: Ord> Heap<T> {

    fn new() -> Heap<T> {
        Heap { elems: Vec::new() }
    }

    fn from(elems: Vec<T>) -> Heap<T> {
        let mut heap = Heap { elems: elems };
        for i in (0..heap.len()/2).rev() {
            heap.max_heapify(i)
        }
        heap
    }

    fn parent(i: usize) -> usize {
        if i > 0 { (i-1)/2 } else { 0 }
    }

    fn left(i: usize) -> usize {
        i*2+1
    }

    fn right(i: usize) -> usize {
        i*2+2
    }

    fn max_heapify(&mut self, i: usize) {
        let (left, right, mut largest) = (Heap::<T>::left(i), Heap::<T>::right(i), i);
        if left < self.len() && self.elems[left] > self.elems[largest] {
            largest = left;
        }
        if right < self.len() && self.elems[right] > self.elems[largest] {
            largest = right;
        }
        if largest != i {
            self.elems.swap(largest, i);
            self.max_heapify(largest);
        }
    }

    fn push(&mut self, v: T) {
        self.elems.push(v);
        let mut i = self.elems.len()-1;
        while i > 0 && self.elems[Heap::<T>::parent(i)] < self.elems[i] {
            self.elems.swap(i, Heap::<T>::parent(i));
            i = Heap::<T>::parent(i);
        }
    }

    fn pop(&mut self) -> Option<T> {
        if self.is_empty() {
            None
        } else {
            let b = self.elems.len()-1;
            self.elems.swap(0, b);
            let v = Some(self.elems.pop().unwrap());
            if !self.is_empty() {
                self.max_heapify(0);
            }
            v
        }
    }

    fn is_empty(&self) -> bool {
        self.elems.is_empty()
    }

    fn len(&self) -> usize {
        self.elems.len()
    }
}

#[cfg(test)]
mod heap_tests {

    use Heap;

    #[test]
    fn test_create() {
        let mut heap: Heap<i32> = Heap::from(vec!(1,3,0,7,8,9,5,6));
        assert_eq!(Some(9), heap.pop());
        assert_eq!(Some(8), heap.pop());
        assert_eq!(Some(7), heap.pop());
        assert_eq!(Some(6), heap.pop());
        assert_eq!(Some(5), heap.pop());
        assert_eq!(Some(3), heap.pop());
        assert_eq!(Some(1), heap.pop());
        assert_eq!(Some(0), heap.pop());
        assert_eq!(None, heap.pop());
    }

    #[test]
    fn test_push() {
        let mut heap: Heap<i32> = Heap::new();
        for &i in [1,3,0,7,8,9,5,6].iter() {
            heap.push(i)
        }
        assert_eq!(Some(9), heap.pop());
        assert_eq!(Some(8), heap.pop());
        assert_eq!(Some(7), heap.pop());
        assert_eq!(Some(6), heap.pop());
        assert_eq!(Some(5), heap.pop());
        assert_eq!(Some(3), heap.pop());
        assert_eq!(Some(1), heap.pop());
        assert_eq!(Some(0), heap.pop());
        assert_eq!(None, heap.pop());
    }

}

fn heap_sort(nums: &mut Vec<usize>) {
    let mut heap: Heap<usize> = Heap::from(nums.clone());
    for i in (0..nums.len()).rev() {
        nums[i] = heap.pop().unwrap();
    }
}

fn main() {

    // 生成随机数组
    let n: usize = 1000000;
    let mut rng = rand::thread_rng();
    let mut nums: Vec<usize> = Vec::new();
    for _i in 0..n {
        nums.push(rng.gen_range(0, n));
    }

    // 注册排序算法
    let algorithms = vec!(
        SortAlgorithm { name: String::from("Heap Sort"), func: heap_sort },
        SortAlgorithm { name: String::from("Bubble Sort"), func: bubble_sort },
        SortAlgorithm { name: String::from("Insert Sort"), func: insert_sort },
        SortAlgorithm { name: String::from("Selection Sort"), func: selection_sort },
        SortAlgorithm { name: String::from("Shell Sort"), func: shell_sort },
        SortAlgorithm { name: String::from("Merge Sort"), func: merge_sort },
        SortAlgorithm { name: String::from("Quick Sort"), func: quick_sort },
        SortAlgorithm { name: String::from("Count Sort"), func: count_sort }
    );

    // 测试排序算法
    for algo in algorithms.iter() {
        let mut temp_nums = nums.clone();
        let start = Instant::now();
        (algo.func)(&mut temp_nums);
        let duration = start.elapsed();
        println!("{}\t{:?}", algo.name, duration);
        assert!(is_sorted(&temp_nums));
    }

}