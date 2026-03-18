def neuro_score(data):
    raw = (
        data['antioxidant'] * 3 +
        data['anti_inflammatory'] * 3 +
        data['mitochondrial_support'] * 2 +
        data['aggregation_modulation'] * 2
    )
    return min(100, round((raw / 40) * 100))


def score_color(score):
    if score >= 70:
        return "green"
    elif score >= 40:
        return "orange"
    else:
        return "red"


def bbb_color(bbb):
    mapping = {
        "High": "green",
        "Medium": "orange",
        "Low": "red",
        "Low-Med": "orange"
    }
    return mapping.get(bbb, "gray")


def disease_color(level):
    mapping = {
        "High": "green",
        "Med": "orange",
        "Low": "red"
    }
    return mapping.get(level, "gray")


def score_label(score):
    if score >= 70:
        return "Strong"
    elif score >= 40:
        return "Moderate"
    else:
        return "Limited"
