package pl.poznan.put.matching;

import pl.poznan.put.pdb.analysis.MoleculeType;
import pl.poznan.put.pdb.analysis.PdbChain;
import pl.poznan.put.pdb.analysis.PdbCompactFragment;
import pl.poznan.put.pdb.analysis.PdbModel;
import pl.poznan.put.pdb.analysis.PdbResidue;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public final class SelectionFactory {
  public static StructureSelection create(final String name, final PdbModel structure) {
    return SelectionFactory.create(name, structure.getChains());
  }

  public static StructureSelection create(final String name, final Iterable<PdbChain> chains) {
    final List<PdbResidue> residues = SelectionFactory.getAllResidues(chains);
    return StructureSelection.divideIntoCompactFragments(name, residues);
  }

  public static StructureSelection create(
      final String name, final PdbModel structure, final SelectionQuery... selectionQueries) {
    final List<PdbCompactFragment> compactFragments = new ArrayList<>(selectionQueries.length);

    for (final SelectionQuery selectionQuery : selectionQueries) {
      final PdbCompactFragment compactFragment = selectionQuery.apply(structure);
      compactFragments.add(compactFragment);
    }

    return new StructureSelection(name, compactFragments);
  }

  private static List<PdbResidue> getAllResidues(final Iterable<PdbChain> chains) {
    final List<PdbResidue> residues = new ArrayList<>();
    for (final PdbChain chain : chains) {
      residues.addAll(SelectionFactory.getAllResidues(chain));
    }
    return residues;
  }

  private static Collection<PdbResidue> getAllResidues(final PdbChain chain) {
    final List<PdbResidue> chainResidues = chain.getResidues();
    final Collection<PdbResidue> residues = new ArrayList<>(chainResidues.size());

    for (final PdbResidue residue : chainResidues) {
      if (residue.getMoleculeType() != MoleculeType.UNKNOWN) {
        residues.add(residue);
      }
    }
    return residues;
  }

  private SelectionFactory() {
    super();
  }
}
