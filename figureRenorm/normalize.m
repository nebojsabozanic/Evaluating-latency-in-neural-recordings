% z-score: renormalize all profiles
function normalized = normalize(profile)
    normalized = (profile - mean(profile))./std(profile);
end