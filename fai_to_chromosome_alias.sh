awk 'BEGIN {OFS="\t"} 
{
  # Declare an array to keep track of used aliases
  if ($1 ~ /^[0-9]+$/) {
    # For non-numeric or unhandled names, use line number
    print NR, "Chr" NR
  } else {
    # Extract the first character
    chrAlias = "Chr" toupper(substr($1, 1, 1))

    # Check if the alias already exists
    if (!(chrAlias in seen)) {
      seen[chrAlias] = 1
      print chrAlias, $1
    } else {
      # Find the next available alias ID
      while ((chrAlias seen[chrAlias]) in seen) {
        seen[chrAlias]++
      }
      aliasWithID = chrAlias seen[chrAlias]
      print aliasWithID, $1
    }
  }
}' input.fai > chromosome_alias.tab