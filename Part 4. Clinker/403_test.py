from pathlib import Path

folder = Path("/Users/anthonycasiday/Documents/Software/plasmid_analyzer_anthony/Part 4. Clinker/output_folder/3kb clinkers from XerC 1kb clusters")
src = folder / "Cluster_1.html"        # read only
dst = folder / "Cluster_1 copy.html"   # output

text = src.read_text()

# 1. Expose Clinker's internal h (locus helpers incl. flip) and d (chart state/update)
anchor = ",y={check:t=>y.checkDomain"
assert text.count(anchor) == 1, f"anchor found {text.count(anchor)} times"
text = text.replace(
    anchor,
    ",_x=(window._h=h,window._d=d),y={check:t=>y.checkDomain",
    1,
)

# 2. After the first plot, flip REV loci the same way the double-click handler does
old = """
} else {
  plot(data)
}
"""
new = """
} else {
  plot(data);
  data.clusters.forEach(cluster => {
    if (cluster.name.includes("_REV_")) {
      cluster.loci.forEach(locus => window._h.flip(locus));
    }
  });
  window._d.update();
}
"""
assert text.count(old) == 1, f"plot block found {text.count(old)} times"
text = text.replace(old, new, 1)

dst.write_text(text)
print("Finished")