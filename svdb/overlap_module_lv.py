from __future__ import absolute_import

def precise_overlap(chrApos_query, chrBpos_query, chrApos_db, chrBpos_db, distance):
    Adist = abs(chrApos_query - chrApos_db)
    Bdist = abs(chrBpos_query - chrBpos_db)
    if max([Adist, Bdist]) <= distance:
        return max([Adist, Bdist]), True
    return False, False


def isSameVariation(chrApos_query, chrBpos_query, chrApos_db, chrBpos_db, ratio, distance):
    if abs(chrApos_query - chrApos_db) <= distance and abs(chrBpos_query - chrBpos_db) <= distance:

        region_start = min([chrApos_db, chrApos_query])
        overlap_start = max([chrApos_db, chrApos_query])

        region_end = max([chrBpos_db, chrBpos_query])
        overlap_end = min([chrBpos_db, chrBpos_query])

        try:
            event_ratio = float(overlap_end - overlap_start + 1) / \
                float(region_end - region_start + 1)
        except Exception:
            event_ratio = 0

        if event_ratio >= ratio:
            return event_ratio, True
        return None, False
    else:
        return None, False


def variant_overlap(chrA, chrB, chrApos_query, chrBpos_query, chrApos_db, chrBpos_db, ratio, distance):
    match = False
    overlap = False
    if chrA == chrB:
        overlap, match = isSameVariation(
            chrApos_query, chrBpos_query, chrApos_db, chrBpos_db, ratio, distance)
    else:
        overlap, match = precise_overlap(
            chrApos_query, chrBpos_query, chrApos_db, chrBpos_db, distance)
    return overlap, match


def weighted_reciprocal_overlap(chrApos_query, chrBpos_query, chrApos_db, chrBpos_db,
                                 distance_weight=0.3, overlap_weight=0.7):

    dist_A = abs(chrApos_query - chrApos_db)
    dist_B = abs(chrBpos_query - chrBpos_db)
    avg_distance = (dist_A + dist_B) / 2

    max_expected_distance = 10000
    normalized_distance = min(avg_distance / max_expected_distance, 1.0)

    region_start = min(chrApos_db, chrApos_query)
    overlap_start = max(chrApos_db, chrApos_query)
    region_end = max(chrBpos_db, chrBpos_query)
    overlap_end = min(chrBpos_db, chrBpos_query)
    
    overlap_length = max(0, overlap_end - overlap_start + 1)
    region_length = region_end - region_start + 1
    overlap_ratio = overlap_length / region_length if region_length > 0 else 0

    similarity_score = (distance_weight * normalized_distance + overlap_weight * (1 - overlap_ratio))
  
    threshold = 0.5
    is_match = similarity_score < threshold

    return similarity_score, is_match


def levenshtein_distance(seq1, seq2):

    if seq1 is None or seq2 is None or len(seq1) == 0 or len(seq2) == 0:
        return (float('inf'), 1.0)
    
    seq1 = str(seq1).upper()
    seq2 = str(seq2).upper()
    
    len1 = len(seq1)
    len2 = len(seq2)
    max_len = max(len1, len2)
    
    # Create DP table
    dp = [[0] * (len2 + 1) for _ in range(len1 + 1)]
    
    # Initialize base cases
    for i in range(len1 + 1):
        dp[i][0] = i
    for j in range(len2 + 1):
        dp[0][j] = j
    
    # Fill DP table
    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            if seq1[i-1] == seq2[j-1]:
                dp[i][j] = dp[i-1][j-1]  # No operation needed
            else:
                dp[i][j] = 1 + min(
                    dp[i-1][j],      # Deletion
                    dp[i][j-1],      # Insertion
                    dp[i-1][j-1]     # Substitution
                )
    
    raw_distance = dp[len1][len2]
    normalized = raw_distance / max_len if max_len > 0 else 0.0
    
    return (raw_distance, normalized)


def compare_insertion_sequences(seq_query, seq_db, max_edit_distance=0.2):
    if seq_query is None or seq_db is None or len(seq_query) == 0 or len(seq_db) == 0:
        return (None, False)
    
    total_dist, norm_dist = levenshtein_distance(seq_query, seq_db)
    
    is_match = norm_dist <= max_edit_distance
    similarity = 1.0 - norm_dist
    
    return (similarity, is_match)


def insertion_overlap_with_sequence(chrApos_query, chrBpos_query, seq_query, 
                                   chrApos_db, chrBpos_db, seq_db, 
                                   distance, max_edit=0.2):

    pos_dist, pos_match = precise_overlap(chrApos_query, chrBpos_query, chrApos_db, chrBpos_db, distance)
    if not pos_match:
        return (None, False)
    
    seq_similarity, seq_match = compare_insertion_sequences(seq_query, seq_db, max_edit)
    
    if seq_similarity is None:
        return (None, False)
    
    overall_match = pos_match and seq_match

    if overall_match:
        if pos_dist is False:
            pos_dist = 0
        combined_similarity = (seq_similarity + (1.0 - pos_dist / distance if distance > 0 else 0)) / 2
        return (combined_similarity, True)
    else: 
        return (None, False)
