import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import heapq

st.title("✈️ Queue Simulation ✈️")

st.markdown(
    """
    This simulation models how $N$ passengers are served at $M$ airport desks.  

    *Params are in the sidebar! (look at the top left corner if using phone)*
    
    - Each **desk** has its own service time distribution, from which the processing time for incoming passengers is sampled.  
    - Each **passenger** has an individual *speed factor* that modifies their service duration.  
    - Service time is sampled from the desk’s distribution and scaled by the passenger’s speed factor:
    
    $t = N(\mu_{desk}, \sigma_{desk}) * N(1, speed_{passenger})$
    """
)

# --- Sidebar controls ---
with st.sidebar:
    coln, colm, cols = st.columns(3)
    with coln:
        num_passengers = st.slider("N", 10, 300, 100, 10)
    with colm:
        num_desks = st.slider("M", 1, 30, 10, 1)
    with cols:
        speed_factor_std = st.slider("Speed σ", 0.0, 0.3, 0.2, 0.1)

    st.markdown("#### Desks service time params")
    col_mu, col_sig = st.columns(2)
    with col_mu:
        mean_range = st.slider("μ [pseudomins]", 5.0, 15.0, (7.0, 12.0))
    with col_sig:
        std_range = st.slider("σ [pseudomins]", 0.1, 0.8, (0.2, 0.5))
    col1, col2 = st.columns(2)
    with col1:
        color1 = st.color_picker("Color №1", "#3CB6EB")
    with col2:
        color2 = st.color_picker("Сolor №2", "#FF4B4B")
    colors = [color1, color2]
# --- Generate distributions ---
desk_time_params = [(np.random.uniform(*mean_range), np.random.uniform(*std_range)) for _ in range(num_desks)]
passenger_speed_factor = np.clip(np.random.normal(1, speed_factor_std, size=num_passengers), 0.1, None)

# --- Simulation using priority queue ---
timeline = [[] for _ in range(num_desks)]
heap = []

# Init desks
for i in range(num_desks):
    heapq.heappush(heap, (0, i))

for passenger in range(num_passengers):
    min_time = heap[0][0]
    free_desks = []

    while heap and heap[0][0] == min_time:
        _, desk = heapq.heappop(heap)
        free_desks.append(desk)

    desk = np.random.choice(free_desks)

    for other in free_desks:
        if other != desk:
            heapq.heappush(heap, (min_time, other))
    
    mean_time, std_time = desk_time_params[desk]
    service_time = round(np.random.normal(mean_time, std_time) * passenger_speed_factor[passenger], 2)
    start_time = min_time
    end_time = start_time + service_time
    timeline[desk].append((start_time, end_time, service_time, passenger))
    heapq.heappush(heap, (end_time, desk))

# --- Total simulation time ---
total_time = max(end for events in timeline for _, end, _, _ in events)
st.write(f"**Waiting time:** {total_time:.2f} pseudomins")

# --- Gantt chart ---
fig, ax = plt.subplots(facecolor='none')
ax.set_facecolor('none')
ax.tick_params(colors='white')
for spine in ax.spines.values():
    spine.set_color('white')
    
ax.set_xlim(0, total_time)
ax.set_ylim(-0.5, num_desks - 0.5)
ax.set_yticks(range(num_desks))
ax.set_xticks(np.arange(0, total_time + 1, 5*(total_time//100+1)))
ax.set_yticklabels([i+1 for i in range(num_desks)])
ax.set_xlabel("Time", color="white", fontsize=13)
ax.set_ylabel("Desk", color="white", fontsize=13)
ax.set_title("Gantt Chart", color='white', fontsize=15)


for desk, events in enumerate(timeline):
    for i, (start, end, _, passenger) in enumerate(events):
        rect = patches.Rectangle(
            (start, desk - 0.4), end - start, 0.8,
            color=colors[i % len(colors)], alpha=0.8
        )
        ax.add_patch(rect)
        ax.text(start + (end - start) / 2, desk, str(passenger + 1),
                ha='center', va='center', color='white', fontsize=7)
fig.tight_layout()
st.pyplot(fig)
