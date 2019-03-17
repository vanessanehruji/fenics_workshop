---
layout: page
title: Tasks
author: Marin Lauber, Vanessa Nehruji
---

## Workshop Tasks

Here is the list of tasks for the workshop:

<ul class="texts">
  {% for item in site.tasks %}
    <li class="text-title">
	<a href="{{ site.baseurl }}{{ item.url }}">
            {{ item.title }}
	</a>
    </li>
    {% endfor %}
</ul>