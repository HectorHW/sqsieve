use std::{
    cmp::Ordering,
    collections::{HashMap, HashSet},
    ops::AddAssign,
};

use itertools::Itertools;

#[derive(Clone, PartialEq, Eq)]
pub struct SparseRow {
    items: Vec<usize>,
}

impl SparseRow {
    fn add_inpace(&mut self, other: &SparseRow) {
        let mut new_buf = vec![];

        let mut i1 = self.items.iter().cloned().peekable();
        let mut i2 = other.items.iter().cloned().peekable();

        while i1.peek().is_some() || i2.peek().is_some() {
            match (i1.peek(), i2.peek()) {
                (Some(&item), None) => {
                    new_buf.push(item);
                    i1.next();
                }

                (None, Some(&item)) => {
                    new_buf.push(item);
                    i2.next();
                }

                (Some(&n1), Some(&n2)) if n1 < n2 => {
                    new_buf.push(n1);
                    i1.next();
                }

                (Some(&n1), Some(&n2)) if n1 > n2 => {
                    new_buf.push(n2);
                    i2.next();
                }

                (Some(_n1), Some(_n2)) => {
                    //terms cancel each other since we are working mod 2
                    i1.next();
                    i2.next();
                }

                _ => unreachable!(),
            }
        }

        self.items = new_buf;
    }

    pub fn contains(&self, item: usize) -> bool {
        self.items.binary_search(&item).is_ok()
    }

    fn is_zero(&self) -> bool {
        self.items.is_empty()
    }

    fn least_term(&self) -> Option<usize> {
        self.items.first().cloned()
    }
}

impl From<&[(usize, usize)]> for SparseRow {
    fn from(value: &[(usize, usize)]) -> Self {
        let mut buf = value
            .iter()
            .filter_map(
                |(term, count)| {
                    if count % 2 == 1 {
                        Some(*term)
                    } else {
                        None
                    }
                },
            )
            .collect_vec();
        buf.sort_unstable();

        SparseRow { items: buf }
    }
}

impl AddAssign<&SparseRow> for SparseRow {
    fn add_assign(&mut self, rhs: &SparseRow) {
        self.add_inpace(rhs)
    }
}

impl PartialOrd for SparseRow {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(match (self.least_term(), other.least_term()) {
            (None, None) => Ordering::Equal,
            (Some(_), None) => Ordering::Less, //we want zeros to go last
            (None, Some(_)) => Ordering::Greater,
            _ => self.items.cmp(&other.items),
        })
    }
}

impl Ord for SparseRow {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

pub struct CongruenceSystem {
    rows: Vec<SparseRow>,
    row_labels: Vec<usize>,
    x_labels: Vec<usize>,
}

pub type SparseCountMap = Vec<Vec<(usize, usize)>>;

impl CongruenceSystem {
    pub fn new(rows: &SparseCountMap, row_labels: Vec<usize>) -> Self {
        assert_eq!(rows.len(), row_labels.len());
        let rows = rows
            .iter()
            .map(|item| SparseRow::from(&item[..]))
            .collect_vec();

        let mut labels = rows
            .iter()
            .flat_map(|row| row.items.iter().cloned())
            .collect::<HashSet<usize>>()
            .into_iter()
            .collect_vec();

        labels.sort_unstable();

        Self {
            rows,
            x_labels: labels,
            row_labels,
        }
    }

    pub fn with_labels(
        rows: &SparseCountMap,
        x_labels: Vec<usize>,
        row_labels: Vec<usize>,
    ) -> Self {
        assert_eq!(rows.len(), row_labels.len());

        let rows = rows
            .iter()
            .map(|item| SparseRow::from(&item[..]))
            .collect_vec();

        Self {
            rows,
            row_labels,
            x_labels,
        }
    }

    fn reorder_descending_slice(slice: &mut [SparseRow]) {
        slice.sort();
    }

    pub fn reorder_descending(&mut self) {
        Self::reorder_descending_slice(&mut self.rows)
    }

    pub fn diagonalize(&mut self) {
        self.reorder_descending();

        for row_number in 0..self.rows.len() {
            if self.rows[row_number].is_zero() {
                break;
            }

            let term = self.rows[row_number].least_term().unwrap();

            for affected_row in (row_number + 1)..self.rows.len() {
                let Some(other_term) = self.rows[affected_row].least_term()  else {
                    continue;
                };
                if other_term == term {
                    //these two rows are always different
                    let (reference, affected) = self.rows.split_at_mut(affected_row);
                    affected[0] += &reference[row_number];
                }
            }

            Self::reorder_descending_slice(&mut self.rows[(row_number + 1)..])
        }
    }

    pub fn print_as_dense(&self) -> String {
        self.rows
            .iter()
            .map(|row| {
                self.x_labels
                    .iter()
                    .map(|&term| if row.contains(term) { "1" } else { "0" })
                    .collect_vec()
                    .join(" ")
            })
            .zip(&self.row_labels)
            .map(|(s, l)| format!("{l: >4}: {s}"))
            .collect_vec()
            .join("\n")
    }

    pub fn transpose(self) -> CongruenceSystem {
        let new_x = self.row_labels;
        let new_row_labels = self.x_labels;

        let old_rows = self.rows;

        let new_rows = new_row_labels
            .iter()
            .map(|&row_label| {
                let items = old_rows
                    .iter()
                    .zip(new_x.iter())
                    .filter_map(|(a, &b)| {
                        if a.contains(row_label) {
                            Some((b, 1usize))
                        } else {
                            None
                        }
                    })
                    .collect_vec();

                SparseRow::from(&items[..])
            })
            .collect_vec();

        CongruenceSystem {
            rows: new_rows,
            row_labels: new_row_labels,
            x_labels: new_x,
        }
    }

    /// fast pivoting algorithm. Matrix is expected to have a column for each smooth number
    pub fn fast_pivot(&mut self) -> Vec<Vec<usize>> {
        assert!(self.x_labels.len() > self.rows.len());
        let mut marking = HashSet::new();
        let mut pivots = HashMap::new();

        for row_number in 0..self.rows.len() {
            let Some(i) = self.rows[row_number].least_term() else {
                continue;
            };
            pivots.insert(row_number, i);
            marking.insert(i);

            for k in 0..self.rows.len() {
                if k == row_number {
                    continue;
                }

                if self.rows[k].contains(i) {
                    let right_side = self.rows[row_number].clone();
                    self.rows[k].add_inpace(&right_side);
                }
            }
        }

        let mut result = vec![];

        for &number in &self.x_labels {
            if marking.contains(&number) {
                continue;
            }

            let mut subresult = vec![number];
            for row_id in 0..self.rows.len() {
                if self.rows[row_id].contains(number) {
                    subresult.push(pivots[&row_id]);
                }
            }
            result.push(subresult);
        }

        result
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Dependency {
    pub variable: usize,
    pub factors: HashSet<usize>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Solution {
    pub vars: HashSet<usize>,
    /// variables that are used to calculate other variables
    pub free_variables: HashSet<usize>,
    /// variables that do not impact solution
    pub lonely_variables: HashSet<usize>,
    /// variables that are always equal to zero
    pub constants: HashSet<usize>,
    pub dependencies: Vec<Dependency>,
}

impl Solution {
    pub fn subsitute(&self, free_vars: &[bool], one_lonelies: bool) -> Vec<bool> {
        assert_eq!(free_vars.len(), self.free_variables.len());
        let mut answer = vec![false; self.vars.len()];

        for (&idx, &value) in self.free_variables.iter().zip(free_vars.iter()) {
            answer[idx] = value;
        }

        if one_lonelies {
            for &var in &self.lonely_variables {
                answer[var] = true;
            }
        }

        //constants are always zero which is true by creation

        for dep in &self.dependencies {
            let sum = dep
                .factors
                .iter()
                .cloned()
                .map(|idx| answer[idx])
                .reduce(|a, b| a != b)
                .unwrap();
            answer[dep.variable] = sum;
        }

        answer
    }
}

pub fn produce_solution(system: &CongruenceSystem) -> Solution {
    let mut constants = HashSet::new();
    let mut const_row_indices = HashSet::new();

    for (i, row) in system.rows.iter().enumerate() {
        let Some(term) = row.least_term() else {
            break;
        };

        if row.items.len() == 1 {
            constants.insert(term);
            const_row_indices.insert(i);
        }
    }

    let dependencies = system.rows.iter().rev().filter(|row| {
        row.items.len() >= 2 && {
            let items = row.items.iter().cloned().collect::<HashSet<usize>>();
            let deps = items.difference(&constants);
            deps.count() > 1
        }
    });

    let mut dependant_vars = HashSet::new();
    let mut free_variables = HashSet::new();
    let mut relations = vec![];

    for dep in dependencies {
        dependant_vars.insert(dep.least_term().unwrap());
        let right_side = dep
            .items
            .iter()
            .cloned()
            .skip(1)
            .filter(|item| !constants.contains(item))
            .collect::<HashSet<_>>();
        free_variables.extend(
            right_side
                .iter()
                .filter(|item| !dependant_vars.contains(item)),
        );
        relations.push(Dependency {
            variable: dep.least_term().unwrap(),
            factors: right_side,
        });
    }

    let participating = dependant_vars
        .iter()
        .chain(free_variables.iter())
        .chain(constants.iter())
        .cloned()
        .collect::<HashSet<_>>();

    let lonely_variables = system
        .x_labels
        .iter()
        .cloned()
        .filter(|var| !participating.contains(var))
        .collect::<HashSet<_>>();

    let vars = system.x_labels.iter().cloned().collect::<HashSet<usize>>();

    Solution {
        vars,
        free_variables,
        lonely_variables,
        constants,
        dependencies: relations,
    }
}

#[cfg(test)]
mod tests {

    use super::{produce_solution, CongruenceSystem, Dependency, Solution, SparseRow};
    use std::collections::HashSet;

    macro_rules! set {
    ( $( $x:expr ),* ) => {  // Match zero or more comma delimited items
        {
            let mut temp_set = HashSet::new();  // Create a mutable HashSet
            $(
                temp_set.insert($x); // Insert each item matched into the HashSet
            )*
            temp_set // Return the populated HashSet
        }
    };
    }

    #[test]
    fn add_sparse_simple() {
        let mut r1 = SparseRow { items: vec![1, 5] };
        let r2 = SparseRow {
            items: vec![2, 4, 6],
        };
        r1.add_inpace(&r2);
        assert_eq!(r1.items, vec![1, 2, 4, 5, 6])
    }

    #[test]
    fn terms_should_cancel() {
        let mut r1 = SparseRow { items: vec![1, 5] };
        let r2 = SparseRow {
            items: vec![2, 5, 6],
        };
        r1.add_inpace(&r2);
        assert_eq!(r1.items, vec![1, 2, 6])
    }

    #[test]
    fn should_reorder() {
        let r1 = SparseRow { items: vec![1, 5] };
        let r2 = SparseRow {
            items: vec![2, 5, 6],
        };

        assert!(r1 < r2)
    }

    #[test]
    fn should_solve() {
        // system:
        // 1 1 0 0 0
        // 0 0 1 0 0
        let mut system = CongruenceSystem::with_labels(
            &vec![vec![(0, 1), (1, 1)], vec![(2, 1)]],
            vec![0usize, 1, 2, 3, 4],
            vec![0usize, 1],
        );

        system.diagonalize();

        assert_eq!(
            produce_solution(&system),
            Solution {
                vars: set![0usize, 1, 2, 3, 4],
                free_variables: set![1usize],
                lonely_variables: set![3usize, 4],
                constants: set![2usize],
                dependencies: vec![Dependency {
                    variable: 0,
                    factors: set![1]
                }]
            }
        )
    }

    trait Unorder {
        fn unorder(self) -> HashSet<usize>;
    }

    impl Unorder for Vec<usize> {
        fn unorder(self) -> HashSet<usize> {
            self.into_iter().collect()
        }
    }

    #[test]
    fn should_solve_with_fast_pivot() {
        let mut system = CongruenceSystem::with_labels(
            &vec![
                vec![(0, 1), (1, 1)],
                vec![(0, 1), (1, 1), (2, 1)],
                vec![(2, 1), (3, 1)],
                vec![(1, 1), (2, 1), (4, 1)],
            ],
            vec![0, 1, 2, 3, 4],
            vec![2, 3, 5, 7],
        );

        println!("{}", system.print_as_dense());

        let result = system.fast_pivot();
        println!("{:?}", result);
        assert_eq!(result.len(), 1);
        assert_eq!(
            result.into_iter().next().unwrap().unorder(),
            vec![0, 1, 4].unorder()
        )
    }
}
