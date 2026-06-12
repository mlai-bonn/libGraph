# TODO

Findings from a code review on 2026-06-12 (analysis only, no code was changed). Ordered by priority.

## Bugs (high priority)

- [ ] **Dangling references in `GraphStruct` data accessors** — `include/GraphDataStructures/GraphBase.h:153` and `:157`: the base-class defaults `GetNodeData(NodeId)` and `GetEdgeData(const EDGE&)` return a `const std::vector<double>&` to a temporary (`return std::vector<double>{0.0};`). Undefined behavior whenever these defaults are hit (e.g. calling `GetEdgeData` on a plain `GraphStruct` in the GED edit-path code). Return by value or a reference to a static empty vector instead.

- [ ] **Missing `break` causes case fall-through** — `include/Algorithms/GED/GEDEvaluation.h:340-347`: the `RandomInsertEdges` case in `get_edit_path` has no `break` and falls through into `RandomDeleteNodes`, so node deletions get inserted (and `missing_node_deletions` cleared) whenever `RandomInsertEdges` is processed. (`DeleteIsolatedNodes` at line 384 also lacks a `break`, currently harmless as the last case, but fragile.)

- [ ] **Use-after-move in `remove_edge`** — `include/Algorithms/GED/GEDEvaluation.h:485-502`: `new_graph` is `std::move`d into `edit_path.Update(...)`, then `new_graph.degree(...)` is queried for isolated-node handling. Reading a moved-from graph is unreliable; query `edit_path.edit_path_graphs.back()` (or check the degree before the move).

- [ ] **Dead "applied twice" guard / potential infinite loop** — `include/Algorithms/GED/GEDEvaluation.h:403-432`: `last_operation` is declared *inside* the `while` loop, so it is always default-constructed; the duplicate-application check can never fire and the trailing `last_operation = next_operation` is dead. If an operation ever fails to be removed by `Update()` (e.g. the early `return` path in `remove_node` at line 555), the loop spins forever. Hoist `last_operation` out of the loop and/or make each iteration pop the front explicitly.

- [ ] **Off-by-one in deleted/inserted node detection** — repeated pattern in `include/Algorithms/GED/GEDEvaluation.h:61, 70, 715-716` and `include/Algorithms/GED/GEDLIBWrapper.h:148-149`: deleted nodes are detected with `x > nodes()` / `x > mapping.size()`, which misclassifies `x == nodes()` (already out of range, since valid ids are `0..nodes()-1`). Should be `>=`. Same question for the "inserted" checks in `add_edge` (`current_i <= std::max(...)` at lines 591-592, should likely be `<`).

- [ ] **`EditPathStrategiesToStringShort` UB on empty result** — `include/Algorithms/GED/GEDStructs.h:121`: `result.back()` is called without checking `result.empty()`; an empty strategy vector (or one containing only strategies that map to `""`, i.e. all `Random*` variants — see `EditPathStrategyToStringShort` default case at line 104) crashes. Also: the `Random*` variants have no short name at all, so folder names silently lose information.

- [ ] **ODR violations: non-inline free functions in header** — `include/Algorithms/GED/GEDLIBWrapper.h:217, 276, 308, 333`: `GEDMethodToString`, `GEDMethodFromString`, `EditCostsToString`, `EditCostsFromString` are non-template, non-`inline` functions defined in a header. Including the header from two translation units breaks the link. Mark them `inline`.

- [ ] **`const_cast` / raw-pointer lifetime hazard in `GEDEvaluation`** — `include/Algorithms/GED/GEDLIBWrapper.h:117-118` `const_cast`s away the constness of `graphs.graphData[...]`, and `GEDEvaluation<T>::graphs` (`include/Algorithms/GED/GEDEvaluation.h:19`) stores raw non-owning pointers into the `GraphData` container. Results silently dangle if the `GraphData` is resized, moved, or destroyed. Consider storing the graph ids only (already present) and resolving graphs on demand, or using `std::reference_wrapper` with documented lifetime requirements.

- [ ] **Broken `executables/test` target** — `executables/test/CMakeLists.txt` lists `../../include/Algorithms/GED/GEDApproximation.h`, which does not exist (the file is now `GEDFunctions.h`/`GED.h`), and compiles both `test.cpp` and `../../examples/graphs/main.cpp` into one target → two `main()` definitions. The target cannot configure/link as-is.

- [ ] **`ComputeGEDResults` directory logic is contradictory** — `include/Algorithms/GED/GEDLIBWrapper.h:170-181`: the function returns early if `results_path` does *not* exist, then immediately calls `create_directory(results_path)` (dead code). Either create the directory when missing or drop the `create_directory` call.

## Correctness / robustness (medium priority)

- [ ] **`CheckResultValidity` duplicate heuristic** — `include/Algorithms/GED/GEDLIBWrapper.h:146-160` and the near-identical copy in `GEDEvaluation.h:708-729` (`CheckResultsValidity`): the duplicate detection collapses *all* deleted-node sentinels into one set entry and compensates with `std::max(0, deleted-1)`. This silently assumes every deleted node maps to the *same* sentinel value. Deduplicate the logic (one function, reused) and make the sentinel explicit (e.g. `ged::GEDGraph::dummy_node()`).

- [ ] **`counter` semantics in `ComputeGEDResults`** — `include/Algorithms/GED/GEDLIBWrapper.h:183-214`: invalid results are saved but not counted, and the (commented-out) time estimate divides by `counter + 1` while iterating `graph_ids` but estimates against `n*(n-1)/2` pairs even when `graph_ids` is a subset. Clarify what the return value means.

- [ ] **Binary serialization is platform-dependent** — `EditOperation::WriteToBinary/ReadFromBinary` (`include/Algorithms/GED/GEDStructs.h:528-542`) and the `.bgf(s)` writers (`GraphBase.h:1840ff`) write raw `size_t`/enum bytes. Files are not portable between platforms (type width, endianness) and there is no magic number/version in the edit-operation stream. At minimum document this; ideally write fixed-width types.

- [ ] **Signed/unsigned mixing** — recurring pattern: `int`/`auto` loop counters compared against `size_t` (`GEDLIBWrapper.h:72, 81`, `GEDEvaluation.h:91, 127`, `GEDStructs.h:214`), `NodeId node = -1` as sentinel in `EditOperation` (`GEDStructs.h:505-506`, `NodeId` is `size_t`), and `edit_path.source_to_current[source_node] = -1` (`GEDEvaluation.h:575`). Works by wraparound but is fragile and triggers warnings; define an explicit `INVALID_NODE` constant.

- [ ] **Missing direct includes** — `GEDLIBWrapper.h` uses `std::filesystem`, `std::chrono`, `std::set`, `std::count_if`, `std::cout` without including the corresponding headers; `GEDStructs.h` uses `std::cerr`/`std::vector`/`std::list`/`std::fstream` but includes only `<iosfwd>`, `<stdexcept>`, `<unordered_set>`. Currently compiles via transitive includes from `GraphBase.h`; add the direct includes.

- [ ] **`GetValidStrategy` complexity** — `include/Algorithms/GED/GEDStructs.h:195-432`: ~240 lines of copy-pasted per-enum bookkeeping with index-juggling `erase` inside a `while` loop. A `std::set`-based dedup plus a small conflict table would shrink this to ~30 lines and remove the off-by-one risk.

- [ ] **Weak hash combine** — `EditOperationHash` (`include/Algorithms/GED/GEDStructs.h:546-556`) combines hashes with shifted XORs; for small integer ids this collides heavily. Use a standard hash-combine (boost::hash_combine pattern).

- [ ] **Error handling via `std::cerr` + return-value-less continuation** — most validation paths print to `std::cerr` and continue (e.g. `get_edit_path` size checks at `GEDEvaluation.h:391-397`). After the recent change to throw in `StringToEditPathStrategy`, consider making error handling consistent (exceptions or a status return) across the GED module.

## Testing / infrastructure

- [ ] **No tests for the GED module** — `Google_tests/TestsMain.cpp` covers graph classes, algorithms, closures, cores, and outerplanar subgraphs, but nothing under `include/Algorithms/GED/`. The edit-path machinery (strategy ordering, node-id remapping in `source_to_current`/`target_to_current`) is exactly the kind of code that needs unit tests; most bugs above would have been caught by one.
- [ ] **No root `CMakeLists.txt`** — the library cannot be consumed via `add_subdirectory`/`find_package`; every consumer hand-rolls `include_directories(../include)`. Add a root CMake with an `INTERFACE` target (`libGraph::libGraph`) and wire the examples/tests/executables as subprojects.
- [ ] **No CI** — no `.github/workflows`; at minimum build the GoogleTest suite and run it on push.
- [ ] **Hard-coded compiler** — `executables/test/CMakeLists.txt` sets `CMAKE_CXX_COMPILER /usr/bin/g++-11`, which breaks on machines without that exact path; pass via toolchain/environment instead.
- [ ] **Build artifacts and IDE folders tracked in the working tree** — `cmake-build-*`, `executables/test/build/`, `Results/`, `.idea/` are untracked clutter; extend `.gitignore`.
- [ ] **`executables/geodesic_core` has no sources** — only `cmake-build-*` and `results/` remain; either restore the source/CMake project or remove the folder.

## Cleanup (low priority)

- [ ] **Leftover commented-out code** — large dead blocks, e.g. `GEDEvaluation.h:434-466` (old edit-path loop) and many commented `std::cout` progress prints in `GEDLIBWrapper.h`/`GEDFunctions.h`. Remove or put behind a verbosity flag.
- [ ] **`_graph` rename artifacts in docs/comments** — comments like "undirected _graph" / "directed data _graph" (e.g. `GraphDirectedFeaturedBase.h:137`, `GraphBase.h:764`) are leftovers of a member rename applied to prose; fix the wording (README was already corrected).
- [ ] **Redundant `result.valid = true` else-branch** — `GEDLIBWrapper.h:199-205`: `CheckResultValidity` already sets `result.valid`; setting it again in the caller's `else` duplicates state logic. Better: have `CheckResultValidity` set `valid` in both directions and drop the caller's `else`.
- [ ] **`unsigned long fixed_operations` counting is inconsistent** — `GEDEvaluation.h:265-330`: deterministic `RelabelEdges`/`RelabelNodes` cases do not increment `fixed_operations` while the four other deterministic cases do, so later `Random*` insert positions may land inside the supposedly fixed relabel prefix. Decide whether relabels are "fixed" and count them consistently.
- [ ] **Old TODOs in core classes** — long-standing markers worth triaging: missing in/out-degree for directed graphs (`GraphDirectedBase.h:66`), "throw exception" placeholders in load paths (several files), `connected_only` parameter of `get_edit_path` unimplemented (`GEDEvaluation.h:197`).
- [ ] **README/format mismatch risk** — the `.bgf(s)` table says "Repeat for every node feature" where edge features are meant (fixed in README); verify the writer in `GraphBase.h:1840-1875` matches the documented layout, and document the `compatibility_format_version` semantics.
