package main

import (
  "os"
  "io/ioutil"
  "sort"
  "time"
  "strings"
  "strconv"
  "fmt"
  "math"
  "math/rand"
  "flag"
  "bits"
)

type MoleculeExactMass struct {
  smiles string
  id string
  ndx int
  exact_mass float64
  seq int
  group int
  nbrs int
  nbr_ids bits.BitField
}

func (c * MoleculeExactMass) Selected () bool {
  return c.group >= 0
}

func (c * MoleculeExactMass) Exact_Mass () float64 {
  return c.exact_mass
}

func (c * MoleculeExactMass) Set_Group (g int) {
  c.group = g
}

func (c * MoleculeExactMass) Id () string {
  return c.id
}

func (c * MoleculeExactMass) Smiles () string {
  return c.smiles
}

func (c * MoleculeExactMass) IsNbr (id int) bool {
  return c.nbr_ids.Bit(id)
}

// A quick way of updating the number of still active nbrs

func (c * MoleculeExactMass) Item_Has_Been_Selected(i int) bool {
  if c.nbr_ids.Bit(i) {
    c.nbr_ids.ClrBit(i)
    c.nbrs--
    return true
  }

  return false
}

func (c * MoleculeExactMass) build (buffer string,
                                  input_separator string,
                                  n int) bool {
  f := strings.Split(strings.TrimRight(buffer, "\n"), input_separator)
  c.group = -1
  c.smiles = f[0]
  c.id = f[1]
  tmp,err := strconv.ParseFloat(f[2], 64)
  if nil != err {
    fmt.Fprintf(os.Stderr, "Invalid exact mass '%s'\n", strings.TrimRight(buffer, "\n"))
    return false
  }
  c.exact_mass = tmp
  c.seq = n

  return true
}
type ID_Mass struct {
  id int
  mass float64
}

func (c * MoleculeExactMass) debug_print(m [] MoleculeExactMass) {
  nbrs := c.nbr_ids.TrueBitsLoHi(0, c.nbr_ids.MaxBit())

  fmt.Fprintf(os.Stderr, "%s %f %d", c.id, c.exact_mass, c.nbrs)
  for _,x := range(nbrs) {
    mi := m[x]
    fmt.Fprintf(os.Stderr, " %f", mi.exact_mass)
    
  }
  fmt.Fprintf(os.Stderr, "\n")
}

func (g * ID_Mass) set (zid int, zmass float64) {
  g.id = zid
  g.mass = zmass
}

type ByMass_IDM [] ID_Mass
func (a ByMass_IDM) Len() int           { return len(a) }
func (a ByMass_IDM) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a ByMass_IDM) Less(i, j int) bool {
  if 0.0 == a[i].mass && 0.0 == a[j].mass {
    return false
  }

  if 0.0 == a[i].mass {
    return false
  }

  if 0.0 == a[j].mass {
    return true
  }

  return a[i].mass < a[j].mass
}

type ByMass [] *MoleculeExactMass
func (a ByMass) Len () int { return len(a);}
func (a ByMass) Swap (i, j int) { a[i],a[j] = a[j],a[i]}
func (a ByMass) Less (i, j int) bool { return a[i].exact_mass < a[j].exact_mass}

type ByNbrs [] *MoleculeExactMass
func (a ByNbrs) Len () int { return len(a);}
func (a ByNbrs) Swap (i, j int) { a[i],a[j] = a[j],a[i]}
func (a ByNbrs) Less (i, j int) bool { 
//fmt.Fprintf(os.Stderr, "Comparing %d %f %d %d with %d %f %d %d\n", i, a[i].exact_mass, a[i].nbrs, a[i].group,
//                                                                   j, a[j].exact_mass, a[j].nbrs, a[j].group)
  if a[i].group >= 0 && a[j].group >= 0 {          // both already placed
    return false
  }
  if a[i].group >= 0 && a[j].group < 0 {
    return false
  }
  if a[i].group < 0 && a[j].group >= 0 {
    return true
  }

  if a[i].nbrs > a[j].nbrs {
    return true
  }

  if a[i].nbrs < a[j].nbrs {
    return false
  }
  return a[i].exact_mass < a[j].exact_mass
}

type Group struct {
  group_size int
  nsel int
  id_mass [] ID_Mass
}

func (g * Group) initialise (group_size int) {
  g.group_size = group_size
  g.nsel = 0
  g.id_mass = make([] ID_Mass, group_size)

  return
}

func (g * Group) debug_print (output * os.File) {
  fmt.Fprintf(output, "Group with %d items ", g.nsel)
  for i := 0; i < g.nsel; i++ {
    fmt.Fprintf(output, " %.4f", g.id_mass[i].mass)
  }
  fmt.Fprintf(output, "\n")
}

func (g * Group) accept (zid int, zmass float64) {

  g.id_mass[g.nsel].set(zid, zmass)
  g.nsel++

  return
}

func (g * Group) smallest_difference (zmass float64) float64 {
  if 0 == g.nsel {
    return 0.0
  }

  rc := math.Abs(g.id_mass[0].mass - zmass)
  for i := 1; i < g.nsel; i++ {
    d := math.Abs(g.id_mass[i].mass - zmass)
    if d < rc {
      rc = d
    }
  }

  return rc
}

func (g * Group) write_results (grp int, m [] MoleculeExactMass, write_smiles bool, output * os.File) {
  var prev_mass float64

  for i := 0; i < g.nsel; i++ {
    ndx := g.id_mass[i].id
    if write_smiles {
      fmt.Fprintf(output, "%s %s %d %d %.4f", m[ndx].Smiles(), m[ndx].Id(), grp, ndx, g.id_mass[i].mass)
    } else {
      fmt.Fprintf(output, "%d %d %s %f", grp, ndx, m[ndx].Id(), g.id_mass[i].mass)
    }

    if m[ndx].exact_mass != g.id_mass[i].mass {
      fmt.Fprintf(os.Stderr, "Mass mismatch %s, ndx %d, in group %f in molecule %f diff %f\n", m[ndx].id, ndx, g.id_mass[i].mass, m[ndx].exact_mass, g.id_mass[i].mass - m[ndx].exact_mass)
    }

    if i > 0 {
      fmt.Fprintf(output, " %.4f\n", m[ndx].exact_mass - prev_mass)
    } else {
      fmt.Fprintf(output, " %d\n", g.nsel)
    }
    prev_mass = m[ndx].exact_mass
  }
}

func (g * Group) number_within_range (m float64, min_diff float64) int {
  rc := 0
  for i := 0; i < g.nsel; i++ {
    if math.Abs(m - g.id_mass[i].mass) < min_diff {
      rc++;
    }
  }

  return rc
}

func (g * Group) sort () {
  sort.Sort(ByMass_IDM(g.id_mass))
}

func (g * Group) ok (m[]MoleculeExactMass) bool {
  for i := 0; i < g.nsel; i++ {
    ndx := g.id_mass[i].id
    if math.Abs(g.id_mass[i].mass - m[ndx].exact_mass) > 0.0001 {
      fmt.Fprintf(os.Stderr, "ok:mass mismatch %s, ndx %d, in group %f in molecule %f diff %f\n", m[ndx].id, ndx, g.id_mass[i].mass, m[ndx].exact_mass, g.id_mass[i].mass - m[ndx].exact_mass)
    }
  }

  return true
}

/*
  We are trying to place an unplaced molecule. Can two groups do some swapping to 
  accommodate this
*/

func (g1 * Group) swap_to_accommodate (m MoleculeExactMass, g2 Group) bool {
  return true
}

func (g * Group) full() bool {
  return g.nsel == g.group_size
}

func (g * Group) size() int {
  return g.group_size
}

func bin_search_by_mass (by_mass []*MoleculeExactMass,
                         lhs int,
                         rhs int,
                         m float64) int {
  n := len(by_mass)

  for i := 0; i < n/8; i++ {
    mid := (lhs + rhs) / 2
    if m < by_mass[mid].exact_mass  {
      rhs = mid
    } else {
      lhs = mid
    }
  }

  return lhs
}

func notify_item_selected(m [] MoleculeExactMass,
                          isel int) {
  nbrs := m[isel].nbr_ids.TrueBitsLoHi(0, m[isel].nbr_ids.MaxBit())

  n := len(nbrs)

  for i := 0; i < n; i++ {
    j := nbrs[i]
    m[j].Item_Has_Been_Selected(isel)
  }

  return
}

// To avoid passing too many arguments to Choose_Group we build a struct to 
// hold most of them

type CG_Args struct {
  group_number int
  min_diff float64
  recompute_nbrs int
  placed int
  randomise_best int
  report int
  number_to_select int
}

func (g * Group) Choose_Group(m       [] MoleculeExactMass,
                              by_nbrs []*MoleculeExactMass,
                              args CG_Args,
                              msdiff [] MSDiff)  int {
  n := len(by_nbrs)

//fmt.Fprintf(os.Stderr, "Begin Choose_Group %d\n", args.group_number)
//for k := 0; k < 20; k ++ {
//  fmt.Fprintf(os.Stderr, " %d %s %f %d grp %d\n", k, by_nbrs[k].id, by_nbrs[k].exact_mass, by_nbrs[k].nbrs, by_nbrs[k].group)
//}
  ndx := 0

  rc := 0

  for ; ndx < n && ndx < len(by_nbrs); ndx++ {
    if by_nbrs[ndx].group >= 0 {       // already selected
      continue
    }

//  if ! g.can_take_main(by_nbrs[ndx].ndx, by_nbrs[ndx].exact_mass, args.min_diff, msdiff) {
    if ! g.can_take_id(m, by_nbrs[ndx].ndx) {
      continue
    }

    by_nbrs[ndx].Set_Group(args.group_number)

    args.placed++

    rc++

//  fmt.Fprintf(os.Stderr, "Group %d got %s had %d nbrs\n", args.group_number, by_nbrs[ndx].id, by_nbrs[ndx].nbrs)

    notify_item_selected(m, by_nbrs[ndx].ndx)

    if g.full() {
      return g.size()
    }

    if 0 == args.placed % args.recompute_nbrs {
      sort.Sort(ByNbrs(by_nbrs))
//    fmt.Fprintf(os.Stderr, "%d sorting. selected item %d %s, mass %f nbrs %d\n", args.placed, ndx, by_nbrs[ndx].id, by_nbrs[ndx].exact_mass, by_nbrs[ndx].nbrs)
//    for k := 0; k < 20; k ++ {
//      fmt.Fprintf(os.Stderr, " %d %s %f %d grp %d\n", k, by_nbrs[k].id, by_nbrs[k].exact_mass, by_nbrs[k].nbrs, by_nbrs[k].group)
//    }
      if args.randomise_best > 0 {
        ianshuffle(by_nbrs, args.randomise_best)
      }
      ndx = -1
    }
  }

  return rc
}

func (g * Group) Place_Via_Distance_From_Leader(m [] MoleculeExactMass,
                                by_nbrs []*MoleculeExactMass,
                                by_mass []*MoleculeExactMass,
                                args CG_Args) int {
  n := len(by_nbrs)

  leader := -1

  for i := 0; i < n; i++ {
    if by_nbrs[i].group >= 0 {
      continue
    }
    g.accept(by_nbrs[i].ndx, by_nbrs[i].exact_mass)
    by_nbrs[i].Set_Group(args.group_number)
//  fmt.Fprintf(os.Stderr, "Leader identified, %d mass %f, put in group %d\n", by_nbrs[i].ndx, by_nbrs[i].exact_mass, args.group_number)

    leader = by_nbrs[i].ndx
    break
  }

  if leader < 0 {         // should not happen
    return 0
  }

  rc := -4

//fmt.Fprintf(os.Stderr, "Leader is %d, group %d\n", leader, args.group_number)

// leader is an index in the by_nbrs array. Find it in the by_mass array

  for i := 0; i < len(by_mass); i++ {
    if leader == by_mass[i].ndx {
      leader = i
//    fmt.Fprintf(os.Stderr, "Found leader at by_mass %d, %f %d\n", i, by_mass[i].exact_mass, by_mass[i].group)
      rc = 1
      break
    }
  }

  if rc < 0 || by_mass[leader].group != args.group_number {
    fmt.Fprintf(os.Stderr, "Huh, leader %d (%f) in group %d, but processing group %d\n", leader, by_mass[leader].exact_mass, by_mass[leader].ndx, args.group_number)
    os.Exit(3)
  }

  rc = 1

  for i := 1; i < len(by_mass) && (args.placed + rc) < args.number_to_select; i++ {
    j := leader + i
    if j < len(by_mass) && by_mass[j].group < 0 && g.can_take_id(m, by_mass[j].ndx) {
      by_mass[j].Set_Group(args.group_number)
      notify_item_selected(m, by_mass[j].ndx)
      rc++
    } 

    j = leader - i
    if j >= 0 && by_mass[j].group < 0 && g.can_take_id(m, by_mass[j].ndx) {
      by_mass[j].Set_Group(args.group_number)
      notify_item_selected(m, by_mass[j].ndx)
      rc++
    }
  }

  args.placed += rc

  return rc;
}

func (g * Group) could_take(m [] MoleculeExactMass, id int, mass float64) (bool, float64) {
  if g.group_size == g.nsel {
    return false,0.0
  }

  if mass != m[id].exact_mass {
    fmt.Fprintf(os.Stderr, "could_take:mass mismatch %f vs %f\n", mass, m[id].exact_mass)
  }

  shortest_distance := math.MaxFloat64
  rc := false
  
  for i := 0; i < g.nsel; i++ {
    j := g.id_mass[i].id

    if m[j].IsNbr(id) {
      return false,0.0
    }

    d := math.Abs(mass - g.id_mass[i].mass)

    if d < shortest_distance {
      shortest_distance = d
      rc = true
    }
  }

  return rc, shortest_distance
}

type PlateSelector interface {
  next_group(g [] Group, m float64, min_diff float64, msdiff [] MSDiff) int
}

/*
  Generic function that takes a function that will return the next group
  that will accept a given mass
*/

type plateselector func(m [] MoleculeExactMass, g [] Group, id int) int

func do_place_item(m       [] MoleculeExactMass,
                   by_nbrs []*MoleculeExactMass,
                   g [] Group,
                   cg_args CG_Args, 
                   selector plateselector) int {
  rc := 0
  next_report := len(by_nbrs) + 1
  if cg_args.report > 0 {
    next_report = cg_args.report
  }

  nmolecules := len(by_nbrs)

  ndx := 0
  for ; ndx < len(by_nbrs) && rc < cg_args.number_to_select; ndx++ {
    if by_nbrs[ndx].group >= 0 {
      continue
    }

//  fmt.Fprintf(os.Stderr, "Trying to place %d %s\n", ndx, by_nbrs[ndx].id)
    best_group := selector(m, g, by_nbrs[ndx].ndx)

    if best_group < 0 {
      continue
    }

//  fmt.Fprintf(os.Stderr, "Best group for %f is %d\n", by_nbrs[ndx].exact_mass, best_group)
//  if g[best_group].can_take_main(by_nbrs[ndx].ndx, by_nbrs[ndx].exact_mass, cg_args.min_diff, msdiff) {
    if g[best_group].can_take_id(m, by_nbrs[ndx].ndx) {
      by_nbrs[ndx].Set_Group(best_group)
      notify_item_selected(m, by_nbrs[ndx].ndx)
    } else {
      fmt.Fprintf(os.Stderr, "Very strange, selected group %d did not accept %f\n", best_group, by_nbrs[ndx].exact_mass)
      g[best_group].debug_print(os.Stderr)
      os.Exit(1)
      continue
    }

    rc++

    if 0 == rc % cg_args.recompute_nbrs {
      by_nbrs = remove_selected_items(by_nbrs)
      sort_and_randomise(by_nbrs, cg_args.randomise_best)
      ndx = -1    // back to start of list
    }

    if rc >= next_report {
      fmt.Fprintf(os.Stderr, "Placed %d of %d molecules %.2f\n", rc, nmolecules, float64(rc)/float64(nmolecules))
      next_report += cg_args.report
    }
  }

  return rc
}

func closest_group(m [] MoleculeExactMass,
                   g [] Group,
                   ndx int) int {
  ng := len(g)
  best_group := -1
  shortest_distance := 0.0
  first_unselected_group := -1
  mass := m[ndx].exact_mass

  for i := 0; i < ng; i++ {
    if 0 == g[i].nsel {
      if first_unselected_group < 0 {
        first_unselected_group = i
      }
      continue
    } 

    could_take,d := g[i].could_take(m, ndx, mass)

    if ! could_take {
      continue
    }

    if d < shortest_distance || best_group < 0 {
      best_group = i
      shortest_distance = d
    }
  }

  if first_unselected_group >= 0 {    // may want to revisit this...
    return first_unselected_group
  }

  return best_group
}

func least_occupied_group(m [] MoleculeExactMass,
                          g [] Group,
                          id int) int {
  ng := len(g)
  best_group := -1
  lowest_occupancy := math.MaxInt32

  for i := 0; i < ng; i++ {
    if 0 == g[i].nsel {
      return i
    }

    if g[i].nsel > lowest_occupancy {
      continue
    }

    could_take,_ := g[i].could_take(m, id, m[id].exact_mass)

    if could_take {
      lowest_occupancy = g[i].nsel
      best_group = i
    }
  }

  return best_group
}

func first_group_to_accept(m [] MoleculeExactMass,
                          g [] Group,
                          id int) int {
  ng := len(g)

  for i := 0; i < ng; i++ {
    if 0 == g[i].nsel {
      return i
    }

    could_take,_ := g[i].could_take(m, id, m[id].exact_mass)

    if could_take {
      return i
    }
  }

  return -1
}

// Example is sodium. We want to avoid pairs that differ by a mass and a tolerance

type MSDiff struct {
  mass float64
  diff float64
}

func (msdiff * MSDiff) ok_difference (diff float64) bool {

//fmt.Fprintf(os.Stderr, "Checking diff %f\n", diff)

  if diff < msdiff.mass - msdiff.diff {
    return true
  } else if diff > msdiff.mass + msdiff.diff {
    return true
  }
//fmt.Fprintf(os.Stderr, "ok_difference rejecting %f\n", diff)
  return false
}

func (msdiff * MSDiff) update_in_range (by_mass []*MoleculeExactMass, 
                                        inc bool) {
  n := len(by_mass)

  e1 := msdiff.mass - msdiff.diff
  e2 := msdiff.mass + msdiff.diff

  for i := 0; i < n; i++ {
    mi := by_mass[i].exact_mass
    ndxi := by_mass[i].ndx
    for j := i + 1; j < n; j++ {
      d := by_mass[j].exact_mass - mi
      if d < e1 {    // have not reached the region yet
        continue
      } else if d > e2 {    // we have passed out of the region
        break
      } else {
        ndxj := by_mass[j].ndx
        if inc {
          by_mass[i].nbrs++
          by_mass[j].nbrs++
          by_mass[i].nbr_ids.SetBit(ndxj)
          by_mass[j].nbr_ids.SetBit(ndxi)
        } else {
          by_mass[i].nbrs--
          by_mass[j].nbrs--
          by_mass[i].nbr_ids.ClrBit(ndxj)
          by_mass[j].nbr_ids.ClrBit(ndxi)
        }
//      fmt.Fprintf(os.Stderr, "From %d %f to %d %f increment nbrs (%.4f)\n", i, mi, j, by_mass[j].exact_mass, by_mass[j].exact_mass - mi)
      }
    }
  }

  return
}

func ok_difference (diff float64, msdiff [] MSDiff) bool {
  for i := 0; i < len(msdiff); i++ {
    if ! msdiff[i].ok_difference(diff) {
      return false
    }
  }

  return true
}

func (g * Group) can_take_id (m [] MoleculeExactMass, id int) bool {
  if g.nsel == g.group_size {
    return false
  }

  for i := 0; i < g.nsel; i++ {     // loop through all items currently in group
    j := g.id_mass[i].id
    if m[j].IsNbr(id) {
      return false
    }
  }

  g.id_mass[g.nsel].id = id
  g.id_mass[g.nsel].mass = m[id].exact_mass
  g.nsel++

  return true
}

func remove_selected_items(by_x [] *MoleculeExactMass) [] * MoleculeExactMass {
  ndx := 0
  n := len(by_x)

  for i := 0; i < n; i++ {
    if by_x[i].group < 0 {
      by_x[ndx] = by_x[i]
      ndx++
    }
  }

  if 0 == ndx {
    return by_x
  }

  return by_x[0:ndx]
}

func determine_nbr_counts_sodium(by_mass []*MoleculeExactMass, 
                          min_diff float64,
                          msdiff [] MSDiff,
                          verbose bool) {
  determine_nbr_counts(by_mass, min_diff, verbose)

  for _,x := range(msdiff) {
    x.update_in_range(by_mass, true)
  }

  return
}

func determine_nbr_counts(by_mass []*MoleculeExactMass, 
                          min_diff float64,
                          verbose bool) {

  nmolecules := len(by_mass)
  for i := 0; i < nmolecules; i++ {

    if by_mass[i].Selected() {
      continue
    }

    ndxi := by_mass[i].ndx

    for j := i+1; j < nmolecules; j++ {
      if by_mass[j].Selected() {
        continue
      }
      if by_mass[j].exact_mass - by_mass[i].exact_mass > min_diff {
        break
      } else {
        ndxj := by_mass[j].ndx
        by_mass[i].nbr_ids.SetBit(ndxj)
        by_mass[j].nbr_ids.SetBit(ndxi)
        by_mass[i].nbrs++
        by_mass[j].nbrs++
      }
    }

  }

  if verbose {
    for i,x := range(by_mass) {
      fmt.Fprintf(os.Stderr, "%s ndx %d mass %f %d nbrs\n", x.id, i, x.Exact_Mass(), x.nbrs)
    }
  }
  return
}

func ianshuffle(by_x []*MoleculeExactMass, randomise int) {
//fmt.Fprintf(os.Stderr, "Request to randomise the first %d items\n", randomise)
  if randomise >= len(by_x) {
    return
  }

  for i := 0; i < randomise; i++ {
    if by_x[i].Selected() {
      randomise = i - 1
      break
    }
  }

//fmt.Fprintf(os.Stderr, "Will randomise %d items\n", randomise)
  for i := 0; i < randomise/2; i++ {
    j := rand.Intn(randomise)
//  fmt.Fprintf(os.Stderr, "Swapping %d and %d\n", i, j)
    if j != i {
      by_x[i],by_x[j] = by_x[j],by_x[i]
    }
  }
}

func sort_and_randomise(by_nbrs []*MoleculeExactMass, randomise int) {
  sort.Sort(ByNbrs(by_nbrs))
  if 0 == randomise {
    return
  }

  ianshuffle(by_nbrs, randomise)

  return
}

func file_size (fname string, print_failure_messages bool) int64 {
  st,err := os.Stat(fname)
  if nil != err {
    if print_failure_messages {
      fmt.Fprintf(os.Stderr,"Missing file '%s'\n",fname)
    }
    return 0
  }
  return st.Size()
}

func dash_s (fname string, print_failure_messages bool) bool {
  st,err := os.Stat(fname)
  if nil != err {
    if print_failure_messages {
      fmt.Fprintf(os.Stderr,"Missing file '%s'\n",fname)
    }
    return false
  }
  if 0 == st.Size() {
    if print_failure_messages {
      fmt.Fprintf(os.Stderr,"Empty file '%s'\n", fname)
    }
    return false
  }

  return true
}

func usage (rc int) {
  fmt.Fprintln(os.Stderr, "Assigns molecules to groups without conflicting exact mass values")
  fmt.Fprintln(os.Stderr, "Input is assumed to be a smiles file 'smiles id exact_mass'")
  fmt.Fprintln(os.Stderr, " -g <n>            size of group (default 10)")
  fmt.Fprintln(os.Stderr, " -ng <n>           number of groups to create")
  fmt.Fprintln(os.Stderr, " -m <diff>         minimum exact mass difference within a group (default 3)")
  fmt.Fprintln(os.Stderr, " -i <char>         input file separator (default space, 'tab' is recognised)")
  fmt.Fprintln(os.Stderr, " -U <fname>        write molecules not placed to <fname>")
  fmt.Fprintln(os.Stderr, " -sort             when results are written, sort by mass within each group")
  fmt.Fprintln(os.Stderr, " -smi              write a smiles file")
  fmt.Fprintln(os.Stderr, " -X mass,diff      extra mass exclusion rules")
  fmt.Fprintln(os.Stderr, " -E                create extra groups with items otherwise unable to be placed")
  fmt.Fprintln(os.Stderr, " -r <n>            recompute the density every <n> selections\n")
  fmt.Fprintln(os.Stderr, " -R <n>            after each density calculation, randomise the first <n> items")
  fmt.Fprintln(os.Stderr, " -p <n>            report progress every <n> items processed")
  fmt.Fprintln(os.Stderr, " -G                do selections on a per group basis, default is per item")
  fmt.Fprintln(os.Stderr, " -L                do selections on a per group basis, but use distance from leader")
  fmt.Fprintln(os.Stderr, " -C                do selection on a per molecule basis, but add to closest group")
  fmt.Fprintln(os.Stderr, " -Z                do selection on a per molecule basis, but add to least occupied group")
  fmt.Fprintln(os.Stderr, " -o <n>            group number offset - add <n> to group numbers printed")
  fmt.Fprintln(os.Stderr, " -s <n>            number of header records to skip")
  fmt.Fprintln(os.Stderr, " -v                verbose output")
  fmt.Fprintln(os.Stderr, "")
  fmt.Fprintln(os.Stderr, " With    -smi output is 'smiles id grp ndx mass'")
  fmt.Fprintln(os.Stderr, " Without -smi output is 'grp ndx id mass'")
  fmt.Fprintln(os.Stderr, "               where ndx is the index of that molecule in the input set")
  os.Exit(rc)
}

func main () {

  var input_separator string
  var unplaced_fname string
  var min_diff float64
  var group_size int
  var ngroups int
  var opt int
  var compress_masses_within_groups bool
  var do_sort bool
  var write_smiles bool
  var recompute_nbrs int
  var xstring string
  var extra_groups_for_unplaced bool
  var randomise_best int
  var report int
  var by_group bool
  var by_distance_from_leader bool
  var put_in_closest_group bool
  var put_in_least_occupied_group bool
  var group_number_offset int
  var header_records_to_skip int
  var verbose bool

  flag.IntVar     (&group_size,         "g",      10,    "group size (default 10)")
  flag.IntVar     (&ngroups,            "ng",     0,     "number of groups to select (default as many as possible)")
  flag.Float64Var (&min_diff,           "m",      3.0,     "minimum mass difference within a group (default 3)")
  flag.StringVar  (&input_separator,    "i",     " ",    "Input separator (default space)")
  flag.StringVar  (&unplaced_fname,     "U",     "",     "write unplaced ids to <fname>")
  flag.StringVar  (&xstring,            "X",     "",     "extra mass constraints <mass,diff>")
  flag.IntVar     (&opt,                "opt",     0,     "after initial placement perform optimisation")
  flag.BoolVar    (&do_sort,            "sort", false,    "during output, sort groups by mass")
  flag.BoolVar    (&write_smiles,       "smi",  false,    "write a smiles file")
  flag.BoolVar    (&extra_groups_for_unplaced, "E", false, "Create extra groups for otherwise unplaced items")
  flag.IntVar     (&recompute_nbrs,     "r",     0,         "recompute density every <n> steps")
  flag.IntVar     (&randomise_best,     "R",     0,         "after each sort, randomise the first <n> items")
  flag.IntVar     (&report,             "p",     0,         "report progress every <n> items processed")
  flag.BoolVar    (&by_group,           "G",     false,     "do selections group at a time")
  flag.BoolVar    (&by_distance_from_leader, "L",false,     "do selections group at a time, based on distance from leader")
  flag.BoolVar    (&put_in_closest_group, "C",   false,     "by molecule, add to closest group")
  flag.BoolVar    (&put_in_least_occupied_group, "Z",   false,     "by molecule, add to least occupied group")
  flag.IntVar     (&group_number_offset, "o",    0,         "group number offset (use -o 1 to start group numbers with 1)")
  flag.IntVar     (&header_records_to_skip, "s", 0,         "number of header records to skip")
  flag.BoolVar    (&verbose,            "v",     false, "verbose output")

  flag.Parse()

  if "tab" == input_separator {
    input_separator = "\t"
  } else if "space" == input_separator {
    input_separator = " "
  }

  msdiff := make([]MSDiff, 0)

  if len(xstring) > 0 {
    f := strings.Split(xstring, ",")

    for i := 0; i < len(f); i += 2 {
      m,err := strconv.ParseFloat(f[i], 64)
      if nil != err {
        fmt.Fprintln(os.Stderr, "Invalid mass '%s'\n", f[i])
        os.Exit(1)
      }
      d,err := strconv.ParseFloat(f[i+1], 64)
      if nil != err {
        fmt.Fprintln(os.Stderr, "Invalid mass diff '%s'\n", f[i+1])
        os.Exit(1)
      }

      var md MSDiff
      md.mass = m
      md.diff = d
//    md := MSDiff(m, d)

      msdiff = append(msdiff, md)
    }
  }

  nfiles := len(flag.Args())
  if 0 == nfiles {
    fmt.Fprintln(os.Stderr, "Must specify file(s) to process")
    usage(1)
  }

  if nfiles > 1 {
    fmt.Fprintln(os.Stderr, "Sorry, only processes one file at a time")
    usage(1)
  }

  fname := flag.Args()[0]

  if ! dash_s(fname, true) {
    fmt.Fprintln(os.Stderr, "Missing or empty input file '%s'", fname)
    os.Exit(1)
  }

  zdata,err := ioutil.ReadFile(fname)

  if nil != err {
    fmt.Fprintf(os.Stderr, "Cannot read data from %s\n", os.Args[1])
    os.Exit(1)
  }

  f := strings.Split(string(zdata), "\n")
  f = f[header_records_to_skip:len(f) - 1];

  nmolecules := len(f)

  if verbose {
    fmt.Fprintf(os.Stderr, "Read %d records from '%s'\n", nmolecules, fname)
  }

  if 0 == recompute_nbrs {
    recompute_nbrs = nmolecules + 1
  }

  m := make([]MoleculeExactMass, nmolecules)

  for i := 0; i < nmolecules; i++ {
    if ! m[i].build(f[i], input_separator, i) {
      fmt.Fprintf(os.Stderr, "Cannot process '%s'\n", f[i])
      os.Exit(1)
    }
    m[i].ndx = i
  }

  by_mass := make([]*MoleculeExactMass, nmolecules)
  by_nbrs := make([]*MoleculeExactMass, nmolecules)
  for i := 0; i < nmolecules; i++ {
    by_mass[i] = &m[i]
    by_nbrs[i] = &m[i]
  }

  sort.Sort(ByMass(by_mass))

  if len(msdiff) > 0 {
    determine_nbr_counts_sodium(by_mass, min_diff, msdiff, verbose)
  } else {
    determine_nbr_counts(by_mass, min_diff, verbose)
  }

  sort.Sort(ByNbrs(by_nbrs))

  if randomise_best > 0 {
    rand.Seed( int64(time.Now().UTC().UnixNano()) + 17983 * int64(os.Getpid()))
//  fmt.Fprintf(os.Stderr, "Time %d pid %d\n", time.Now().UTC().UnixNano(), os.Getpid())
    ianshuffle(by_nbrs, randomise_best)
//  for i := 0; i < randomise_best; i++ {        // do our own shuffling because ianshuffle does not move unselected items
//    j := rand.Intn(randomise_best)
//    if j != i {
//      by_nbrs[i],by_nbrs[j] = by_nbrs[j],by_nbrs[i]
//    }
//  }
  }

  if verbose {
    for _,m := range(by_nbrs) {
      fmt.Fprintf(os.Stderr, "%s mass %f has %d nbrs\n", m.id, m.exact_mass, m.nbrs)
    }
  }

  number_to_select := nmolecules

  if 0 == ngroups {
    ngroups = nmolecules / group_size
    if 0 == ngroups {
      fmt.Fprintf(os.Stderr, "From %d molecules cannot create groups of size %d\n", nmolecules, group_size)
      os.Exit(1)
    }
    if ngroups * group_size != nmolecules {
      ngroups++
    }

    if verbose {
      fmt.Fprintf(os.Stderr, "Placing %d items into groups of size %d requires %d groups\n", nmolecules, group_size, ngroups)
    }
  }

  g := make([] Group, ngroups)

  for i := 0; i < ngroups; i++ {
    g[i].initialise(group_size)
  }

// Check to see if there are any impossible combinations

  for i,x := range(by_mass) {
    in_range := 0
    mstop := x.exact_mass + min_diff

    for j := i + 1; j < nmolecules; j++ {
      if by_mass[j].exact_mass > mstop {
        break
      } else {
        in_range++
      }
    }

    if in_range > ngroups {
      fmt.Fprintf(os.Stderr, "Items need %d groups, begin %s\n", in_range, x.id)
    }
  }

  isolated := 0
  for i := 1; i < nmolecules; i++ {
    if by_mass[i].exact_mass - by_mass[i-1].exact_mass > min_diff {
      isolated++
    }
  }

  if verbose {
    fmt.Fprintf(os.Stderr, "%d of %d molecules are isolated within %f\n", isolated, nmolecules, min_diff)
  }

  sort.Sort(ByNbrs(by_nbrs))
  if verbose {
    for _,x := range(by_nbrs) {
      x.debug_print(m)
    }
  }

  placed := 0

  var cg_args CG_Args
  cg_args.min_diff = min_diff
  cg_args.recompute_nbrs = recompute_nbrs
  cg_args.placed = 0
  cg_args.randomise_best = randomise_best
  cg_args.report = report
  cg_args.number_to_select = number_to_select

// The compressing masses within groups idea never really worked, and this code
// does not work. The problem is that we also have the objective
// of placing the first molecules in the first groups

  if compress_masses_within_groups {
    for i := 1; i < nmolecules; i++ {
      closest_match := -1
      min_dist := 9000.0
      for j := 0; j < ngroups; j++ {
        d := g[j].smallest_difference(m[i].Exact_Mass())

        if d > min_diff && d < min_dist {
          closest_match = i
          min_dist = d
        }
      }
      if closest_match >= 0 {
        g[closest_match].accept(i, m[i].Exact_Mass())
        placed++
      }
    }
  } else if by_group {

    for i:= 0; i < ngroups && placed < number_to_select; i++ {
      cg_args.group_number = i
      placed += g[i].Choose_Group(m, by_nbrs, cg_args, msdiff)
      by_nbrs = remove_selected_items(by_nbrs)
      if 0 == i % recompute_nbrs {
        sort_and_randomise(by_nbrs, randomise_best)
      }
    }
  } else if by_distance_from_leader {
    for i := 0; i < ngroups; i++ {
      cg_args.group_number = i
      placed += g[i].Place_Via_Distance_From_Leader(m, by_nbrs, by_mass, cg_args)
      by_nbrs = remove_selected_items(by_nbrs)
      if 0 == i % recompute_nbrs {
        sort_and_randomise(by_nbrs, randomise_best)
      }
    }
  } else if put_in_closest_group {
    placed = do_place_item(m, by_nbrs, g, cg_args, closest_group)
  } else if put_in_least_occupied_group {
    placed = do_place_item(m, by_nbrs, g, cg_args, least_occupied_group)
  } else {
    placed = do_place_item(m, by_nbrs, g, cg_args, first_group_to_accept)
//  placed = do_place_by_first_group_to_accept(m, by_nbrs, g, cg_args, msdiff)
  }

  if verbose {
    fmt.Fprintf(os.Stderr, "Placed %d of %d molecules %.2f%s\n", placed, nmolecules, (float64(placed) / float64(nmolecules) * 100.0), "%")
  }

  if opt > 0 {
    placed_during_opt := 0

    for i := 0; i < nmolecules; i++ {
      if m[i].Selected() {
        continue
      }

      swapped := false
      for j := 0; j < ngroups; j++ {
        for k := j + 1; k < ngroups; k++ {
          if g[j].swap_to_accommodate(m[i], g[k]) {
            swapped = true
            break
          }
        }
        if swapped {
          break
        }
      }
      if swapped {
        placed_during_opt ++
      }
    }

    placed += placed_during_opt

    if verbose {
      fmt.Fprintf(os.Stderr, "During optimisation placed an extra %d molecules, tot %d\n", placed_during_opt, placed)
    }
  }

  if ! write_smiles {
    fmt.Fprintf(os.Stdout, "GRP NDX ID Mass Diff\n")
  }

  for i := 0; i < ngroups; i++ {
    if do_sort {
      g[i].sort()
    }
//  fmt.Fprintf(os.Stderr, "Writing group %d\n", i)
    g[i].write_results (i + group_number_offset, m, write_smiles, os.Stdout)
  }

  if verbose {
    for i:= 0; i < ngroups; i++ {
      fmt.Fprintf(os.Stderr, "%d ", i)
      g[i].debug_print(os.Stderr)
    }
  }

  if verbose && placed < nmolecules {
    fmt.Fprintf(os.Stderr, "Only placed %d of %d items (%d unplaced). Built %d groups\n", placed, nmolecules, nmolecules - placed, ngroups)
  }

//fmt.Fprintf(os.Stderr, "Extra %d, placed %d number_to_select %d\n", extra_groups_for_unplaced, placed, number_to_select)
  if extra_groups_for_unplaced && placed < number_to_select {
    neg := number_to_select - placed    // number of extra groups
    eg := make([]Group, neg, neg)
    for i := 0; i < neg; i++ {
      eg[i].initialise(group_size)
    }

    sort.Sort(ByNbrs(by_nbrs))
    if randomise_best > 0 {
      ianshuffle(by_nbrs, randomise_best)
    }

    placed_here := 0

    for i := 0; i < len(by_nbrs) && placed < number_to_select; i++ {
      if by_nbrs[i].Selected() {
        continue
      }
      iplaced := false
      for j := 0; j < neg; j++ {

        if eg[j].can_take_id(m, by_nbrs[i].ndx) {
          fmt.Fprintf(os.Stderr, "Extra group %d took %f\n", j, by_nbrs[i].exact_mass)
          by_nbrs[i].Set_Group(ngroups + j)
          notify_item_selected(m, by_nbrs[i].ndx)
          iplaced = true
          placed_here++
          break
        }
      }
      if iplaced {
        placed++
      }
    }

    extra_groups_used := 0

    for i,g := range(eg) {
//    fmt.Fprintf(os.Stderr, "Writing extra group %d, contains %d items\n", i, g.nsel)
      if 0 == g.nsel {
        break
      }
      extra_groups_used++;

      g.sort()
      g.write_results (ngroups + i + group_number_offset, m, write_smiles, os.Stdout)
    }

    if verbose {
      fmt.Fprintf(os.Stderr, "Placed %d unplaced items in %d extra groups\n", placed_here, extra_groups_used)
    }
  }

  if len(unplaced_fname) > 0 {
    outp,err := os.Create(unplaced_fname)
    if nil != err {
      fmt.Fprintf(os.Stderr, "Cannot open '%s', %v\n", unplaced_fname, err)
      os.Exit(1)
    }
    fmt.Fprintf(outp, "ID NDX\n")

    for i := 0; i < nmolecules; i++ {
      if m[i].Selected() {
        continue
      }
      fmt.Fprintf(outp, "%s %.4f %d\n", m[i].Id(), m[i].Exact_Mass(), i)
    }
  }

  os.Exit(0)
}
